// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // int32_t, uint32_t
#include <cmath> // std::exp, ::sqrt, ::abs
#include <map> // std::map<Key,T>

#include "parallel_potential.hxx"

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "real_space.hxx" // ::grid_t
#include "display_units.h" // Ang, _Ang
#include "geometry_input.hxx" // ::init_geometry_and_grid
#include "mpi_parallel.hxx" // ::init, ::finalize, MPI_COMM_WORLD, ::barrier
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "data_view.hxx" // view2D<>, view3D<>
#include "global_coordinates.hxx" // ::get
#include "inline_math.hxx" // pow2, set, add_product
#include "exchange_correlation.hxx" // ::LDA_kernel
#include "simple_stats.hxx" // ::Stats<>
#include "constants.hxx" // ::pi
#include "recorded_warnings.hxx" // warn
#include "data_list.hxx" // data_list<T>
#include "boundary_condition.hxx" // ::periodic_images
#include "chemical_symbol.hxx" // ::get
#include "simple_timer.hxx" // SimpleTimer
#include "hermite_polynomial.hxx" // Gauss_Hermite_polynomials
#include "sho_tools.hxx" // ::nSHO, ::n1HO
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform
#include "sho_projection.hxx" // ::denormalize_electrostatics, ::renormalize_electrostatics
#include "print_tools.hxx" // printf_vector
#include "energy_contribution.hxx" // ::show, ::TOTAL, ::KINETIC, ::ELECTROSTATIC, ...
#include "energy_contour.hxx" // ::Integrator

#ifdef    HAS_SINGLE_ATOM
    #include "single_atom.hxx" // ::atom_update
#else  // HAS_SINGLE_ATOM
#ifdef    HAS_LIVE_ATOM
    // use the libliveatom library C-interface
    extern "C" {
       #include "single_atom.h" // live_atom_update_
                                // live_atom_is_a_dynamic_library_
                                // live_atom_init_env_
    } // extern "C"
#endif // HAS_LIVE_ATOM
#endif // HAS_SINGLE_ATOM

namespace parallel_potential {

    status_t live_atom_update(
        char const *const what    // selector string (input)
      , int32_t  const natoms     // number of atoms (input)
      , double  *const dp=nullptr // quantities (input/output) double  dp[natoms]
      , int32_t *const ip=nullptr // quantities (input/output) integer ip[natoms]
      , float   *const fp=nullptr // quantities (input)        float   fp[natoms or less]
      , double  *const *const dpp=nullptr // quantities (input/output) double* dpp[natoms]
    ) {
        int32_t stat(0); // status variable for

        static int use{-1};
        if (-1 == use) {
            use = control::get("use.live.atom", 1.);
#ifndef   HAS_SINGLE_ATOM
#ifdef    HAS_LIVE_ATOM
            int32_t is_dynamic{0}; live_atom_is_a_dynamic_library_(&is_dynamic);
            int const echo = (0 == mpi_parallel::rank());
            if (is_dynamic) {
                auto const *const control_file = control::get("control.file", "");
                if (echo > 0) std::printf("# libliveatom.so is linked as dynamic library\n"
                    "# read single_atom.*-controls from +control.file=%s\n", control_file);
                live_atom_init_env_(control_file, &stat);
                if (0 != stat) {
                    warn("+control.file=%s for libliveatom.so, live_atom_init_env_ returned %i", control_file, int(stat));
                }
            } else {
                // We do not need the control file in the case of a static library 
                //       as we share the control.o and recorded_warnings.o objects
                if (echo > 0) std::printf("# libliveatom.a is linked as static library\n");
            } // is_dynamic
#endif // HAS_LIVE_ATOM
#endif // HAS_SINGLE_ATOM
        } // needs init
        if (0 == use) {
            // warn("single_atom::atom_update deactivated", 0); 
            return 0;
        } // 0 == use

#ifdef    HAS_SINGLE_ATOM
        stat = single_atom::atom_update(what, natoms, dp, ip, fp, dpp);
#else  // HAS_SINGLE_ATOM
#ifdef    HAS_LIVE_ATOM
        live_atom_update_(what, &natoms, dp, ip, fp, dpp, &stat);
#else  // HAS_LIVE_ATOM
        static bool warned = false;
        if (!warned) {
            warn("compiled with neither -DHAS_SINGLE_ATOM nor -DHAS_LIVE_ATOM for %d atoms", natoms);
            warned = true;
        } // launch warning only once
#endif // HAS_LIVE_ATOM
#endif // HAS_SINGLE_ATOM

        if (stat) { warn("single_atom::atom_update(%s, natoms=%d, ...) returned status= %i", what, natoms, int(stat)); }
        return stat;
    } // live_atom_update

    inline double rho_Thomas_Fermi( // Thomas-Fermi model density estimate
        double const Ekin // input: local kinetic energy, i.e. E - V(r)
      , double & derivative // result: d/dE rho
    ) {
        if (Ekin <= 0) { derivative = 0; return 0; }
        double constexpr by_pi2 = 1./pow2(constants::pi);
        derivative = by_pi2 * std::sqrt(2*Ekin);
        double const rho = derivative * Ekin * (2./3.);
        return rho;
    } // rho_Thomas_Fermi


    status_t new_density_Thomas_Fermi(
          double rho_new[] // resulting density
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential
        , size_t const nall // number of all grid points
        , MPI_Comm const comm=MPI_COMM_WORLD // grid communicator
        , double const n_electrons=1 // required total number of electrons 
        , double const dV=1 // grid volume element
        , int const echo=0 // log level
        , int const maxiter=99 // maximum number of iterations
        , float const threshold=1e-8f // convergence criterion
    ) {
        double mu = Fermi_level;
        { // scope: correct Fermi level
            simple_stats::Stats<double> s(0);
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                s.add(Vtot[zyx]);
            } // zyx
            mpi_parallel::allreduce(s, comm);
            if (mu <= s.min()) { mu = 0.5*(s.max() + s.min()); }
            mu = std::max(s.min(), mu);
        } // scope

        if (echo > 3) std::printf("# Thomas-Fermi model, starting with chemical potential mu= %g %s\n", mu*eV, _eV);
        double ne{0}, de{0}; // ne: preliminary number of electrons, de: derivative w.r.t. the Fermi level
        int iterations{0};
        bool converged{false};
        do {
            de = 0; ne = 0;
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                double der;
                double const rho = rho_Thomas_Fermi(mu - Vtot[zyx], der);
                ne += rho;
                de += der;
            } // zyx
            ne = mpi_parallel::sum(ne, comm) * dV;
            de = mpi_parallel::sum(de, comm) * dV;
            auto const residual = std::abs(ne - n_electrons);
            if (echo > 5) std::printf("# Thomas-Fermi iteration #%i mu= %g %s --> %g electrons, derivative= %g, residual= %.1e\n",
                                                                iterations, mu*eV, _eV, ne, de, residual);
            converged = (residual < threshold);
            mu += (n_electrons - ne)/std::max(0.1, de); // new mu by Newton step
            ++iterations;
        } while (!converged && iterations < maxiter);

        if (converged) {
            if (echo > 2) std::printf("# Thomas-Fermi converged in %d iterations, mu= %g %s --> %g electrons, derivative= %g\n",
                                                                iterations, mu*eV, _eV, ne, de);
            for (size_t zyx = 0; zyx < nall; ++zyx) {
                double der;
                rho_new[zyx] = rho_Thomas_Fermi(mu - Vtot[zyx], der);
            } // zyx
            Fermi_level = mu;
            return 0;
        } else { // converged
            warn("Thomas-Fermi model did not converge in %d iterations", iterations);
            return iterations;
        } // converged
    } // new_density_Thomas_Fermi


    template<typename real_t>
    void block_average(
          real_t v4[4*4*4] // result
        , real_t const v8[8*8*8] // input
        , simple_stats::Stats<double> *s=nullptr
    ) {
        for (int iz{0}; iz < 4; ++iz) {
        for (int iy{0}; iy < 4; ++iy) {
        for (int ix{0}; ix < 4; ++ix) {
            auto const sum =  v8[((iz*2    )*8 + iy*2    )*8 + ix*2    ]
                            + v8[((iz*2    )*8 + iy*2    )*8 + ix*2 + 1]
                            + v8[((iz*2    )*8 + iy*2 + 1)*8 + ix*2    ]
                            + v8[((iz*2    )*8 + iy*2 + 1)*8 + ix*2 + 1]
                            + v8[((iz*2 + 1)*8 + iy*2    )*8 + ix*2    ]
                            + v8[((iz*2 + 1)*8 + iy*2    )*8 + ix*2 + 1]
                            + v8[((iz*2 + 1)*8 + iy*2 + 1)*8 + ix*2    ]
                            + v8[((iz*2 + 1)*8 + iy*2 + 1)*8 + ix*2 + 1];
            auto const average = sum*0.125; // divide by 2*2*2
            v4[(iz*4 + iy)*4 + ix] = average;
            if (s) s->add(average);
        }}} // ix iy iz
    } // block_average


    double integrate_r2grid(
          double const r2coeff[]
        , int32_t const nr2=4096
        , float const ar2=16.f
    ) {
        double rho{0};
        auto const by_ar2 = 1./ar2;
        for (int ir2{0}; ir2 < nr2; ++ir2) {
            auto const r2 = ir2*by_ar2, r = std::sqrt(r2);
            rho += r2coeff[ir2]*r*.5*by_ar2;
        } // ir2
        return rho;
    } // integrate_r2grid


    template <typename real_t>
    double print_stats( // similar to print_stat in print_tools.hxx
          real_t const values[] // input values
        , size_t const all // how many local values
        , MPI_Comm const comm // MPI communicator for reduction
        , bool const echo // only true for the master rank
        , double const dV=1 // volume element for the integration
        , char const *prefix="" // leading printf messages
        , double const unit=1 // unit conversion factor
        , char const *_unit="" // unit indicator
    ) {
        simple_stats::Stats<double> s(0);
        for (size_t i{0}; i < all; ++i) {
            s.add(values[i]);
        } // i
        mpi_parallel::allreduce(s, comm);
        if (echo) {
         // std::printf("%s grid stats min %g max %g avg %g", prefix, s.min()*unit, s.max()*unit, s.mean()*unit);
            std::printf("%s grid stats [%g, %g +/- %g, %g]", prefix, s.min()*unit, s.mean()*unit, s.dev()*unit, s.max()*unit);
            if (dV > 0) std::printf(" %g electrons", s.sum()*dV*unit);
            std::printf(" %s\n", _unit);
        } // echo
        return s.sum()*dV;
    } // print_stats


    class atom_image_t {
    public:
        atom_image_t() {}
        atom_image_t(double const pos[3], double const Z, int32_t const atom_id=-1, int8_t const shifts[3]=nullptr) {
            set(pos_, 3, pos);
            atom_id_ = atom_id;
            set(shifts_, 3, int8_t(0)); if (shifts) set(shifts_, 3, shifts);
            iZ_ = chemical_symbol::get(Z);
        } // constructor
    public:
        double pos_[3]; // atom image position
        int32_t atom_id_; // atom index, first global, later local
        int8_t shifts_[3]; // periodic image shifts in [-127, 127]
        int8_t iZ_; // atomic number hint, the true atomic number can be retrieved from the full list
    }; // class atom_image_t, 32 Byte


    std::vector<atom_image_t> get_neighborhood( // determine which atoms are relevant for my blocks
          view2D<double> & block_coords // [n_blocks][4], will be filled with values
        , int32_t const n_all_atoms
        , view2D<double> const & xyzZ_all // [n_all_atoms][4] all atomic coordinates and atomic numbers
        , parallel_poisson::parallel_grid_t const & pg
        , real_space::grid_t const & g
        , double const grid_center[3]
        , int const echo=0
        , float const r_cut=16.f // truncation radius of atom-cube interactions, largest radius for smooth core densities, vbar, 9*sigma_cmp
    ) {
        if (echo > 2) std::printf("# search periodic images of atoms in a neighborhood of radius %g %s\n", r_cut*Ang, _Ang);

        std::vector<atom_image_t> atom_images(0); // list of candidates of atomic images relevant for this MPI rank

        auto const n_blocks = pg.n_local();
        if (n_blocks < 1) return atom_images;

        // determine the sphere around an 8x8x8 cube
        auto const h = g.grid_spacings();
        auto const r_circum = 4*std::sqrt(pow2(h[0]) + pow2(h[1]) + pow2(h[2]));

        assert(block_coords.stride() >= 4);

        // determine the center of weight of all local 8x8x8 cubes
        double cow[] = {0, 0, 0};
        auto const local_ids = pg.local_ids();
        for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
            uint32_t ixyz[3]; global_coordinates::get(ixyz, local_ids[ilb]);
            auto *const block_pos = block_coords[ilb]; 
            for (int d{0}; d < 3; ++d) {
                auto const p = h[d]*(ixyz[d]*8. + 4.) - grid_center[d];
                block_pos[d] = p;
                cow[d] += p;
            } // d
        } // ilb
        assert(n_blocks > 0);
        scale(cow, 3, 1./n_blocks);

        // determine the largest distance of a cube from the center of mass
        double max_dist2{0};
        for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
            auto const *const block_pos = block_coords[ilb];
            double dist2{0};
            for (int d{0}; d < 3; ++d) {
                dist2 += pow2(block_pos[d] - cow[d]);
            } // d
            block_coords[ilb][3] = dist2; // store distance^2 from center of mass
            max_dist2 = std::max(max_dist2, dist2);
        } // ilb
        auto const max_dist = std::sqrt(max_dist2);

        auto const me = mpi_parallel::rank(pg.comm());
        if (echo > 2) std::printf("# rank#%i all grid points within %g + %g = %g %s around %g %g %g %s\n",
            me, max_dist*Ang, r_circum*Ang, (max_dist + r_circum)*Ang, _Ang, cow[0]*Ang, cow[1]*Ang, cow[2]*Ang, _Ang);



        { // scope: construct periodic images, gather candidates
            auto const r_search = max_dist + r_circum + r_cut, r2search = pow2(r_search);

            view2D<double> periodic_images;
            view2D<int8_t> periodic_shift;
            auto const n_periodic_images = boundary_condition::periodic_images(periodic_images,
                           g.cell, g.boundary_conditions(), r_search, echo/2, &periodic_shift);

            // this should be the only loop over all atoms in the entire code
            // keep it free of I/O and only minimum number of operations
            {   SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances of all atoms", echo);
                for (int32_t global_atom_id{0}; global_atom_id < n_all_atoms; ++global_atom_id) {
                    auto const *const pos = xyzZ_all[global_atom_id];
                    for (int ip{0}; ip < n_periodic_images; ++ip) {
                        auto const *const img = periodic_images[ip];
                        auto const dist2 = pow2(pos[0] + img[0] - cow[0])
                                         + pow2(pos[1] + img[1] - cow[1])
                                         + pow2(pos[2] + img[2] - cow[2]);
                        if (dist2 < r2search) {
                            double const img_pos[] = {pos[0] + img[0], pos[1] + img[1], pos[2] + img[2]};
                            double const Z = pos[3]; // pos points into xyzZ array
                            atom_images.push_back(atom_image_t(img_pos, Z, global_atom_id, periodic_shift[ip]));
                        } // inside search radius
                    } // ip
                } // global_atom_id
            } // timer
            if (echo > 2) std::printf("# rank#%i finds %ld atom images inside a %g %s search radius\n",
                                              me, atom_images.size(), r_search*Ang, _Ang);
        } // scope


        // now refine the search checking the proximity
        auto const r_close = r_cut + r_circum, r2close = pow2(r_close);

     // view2D<float> atom_block_distance2(atom_images.size(), n_blocks, 9e9);
        std::vector<uint32_t> n_close_blocks(atom_images.size(), 0),
                              n_close_atoms(n_blocks, 0);
        size_t distances_close{0}; // how many atom images are close to the domain

        {   SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances", echo);
            for (size_t ja{0}; ja < atom_images.size(); ++ja) { // OMP PARALLEL
                auto const *const atom_pos = atom_images[ja].pos_;
                uint32_t n_blocks_close{0};
                for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
                    auto const *const block_pos = block_coords[ilb];

                    double dist2{0};
                    for (int d{0}; d < 3; ++d) {
                        dist2 += pow2(block_pos[d] - atom_pos[d]);
                    } // d
                // atom_block_distance2(ja,ilb) = dist2; // store
                    if (dist2 < r2close) {
                        ++n_close_atoms[ilb];
                        ++n_blocks_close;
                    }
                } // ilb
                n_close_blocks[ja] = n_blocks_close;
                distances_close   += n_blocks_close;
            } // ja
        } // timer
        auto const distances_checked = n_blocks*atom_images.size();
        if (echo > 2) std::printf("# rank#%i finds %.6f M of %.6f M atom-block distances inside a %g %s search radius\n",
                                          me, distances_close*1e-6, distances_checked*1e-6, r_close*Ang, _Ang);

        { // scope: check if any blocks are far from atoms
            uint32_t n_far_blocks{0};
            for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
                n_far_blocks += (0 == n_close_atoms[ilb]);
            } // ilb
            if (n_far_blocks && echo > 1) std::printf("# rank#%i has %d blocks far from any atom image\n", me, n_far_blocks);
        } // scope

        // check if any atom images can be removed from the list
        uint32_t n_relevance{0};
        std::vector<int32_t> new_index(atom_images.size(), -1);
        for (size_t ja{0}; ja < atom_images.size(); ++ja) {
            if (n_close_blocks[ja] > 0) {
                new_index[n_relevance] = ja;
                ++n_relevance;
            }
        } // ja
        auto const n_relevant_atoms = n_relevance;
        auto const shrink = (n_relevant_atoms < atom_images.size());
        if (echo > 2) std::printf("# rank#%i %s %.3f k of %.3f k atom images are relevant\n",
                    me, shrink?"only":"all", n_relevant_atoms*1e-3, atom_images.size()*1e-3);

        if (shrink) {
            SimpleTimer timer(strip_path(__FILE__), __LINE__, "removing irrelevant atom images", echo);
            // there are irrelevant atoms, delete them from the list
            std::vector<uint32_t>  n_close_blocks_(n_relevant_atoms);
            std::vector<atom_image_t> atom_images_(n_relevant_atoms);
         // view2D<float> atom_block_distance2_(n_blocks, n_relevant_atoms);
            for (uint32_t ka{0}; ka < n_relevant_atoms; ++ka) {
                auto const ja = new_index[ka];
                assert(ja >= 0);
                atom_images_[ka]    = atom_images[ja];
                n_close_blocks_[ka] = n_close_blocks[ja];
             // set(atom_block_distance2_[ka], n_blocks, atom_block_distance2[ja]);
            } // ka
            n_close_blocks = n_close_blocks_; // overwrite
            atom_images    = atom_images_;    // overwrite
            // atom_block_distance2 = atom_block_distance2_; // this does not compile since the move operators are deleted on purpose
         // atom_block_distance2 = view2D<float>(n_relevant_atoms, n_blocks);
         // for (uint32_t ka{0}; ka < n_relevant_atoms; ++ka) {
         //     set(atom_block_distance2[ka], n_blocks, atom_block_distance2_[ka]); // copy
         // } // ka
        } // delete irrelevant atoms
        new_index.clear(); // not needed any longer
        assert(atom_images.size() == n_relevant_atoms);
        auto const distances_stored = n_blocks*n_relevant_atoms;
     // if (echo > 2) std::printf("# rank#%i %d x %.3f k = %.6f M atom-block distances stored, %.3f MByte\n",
     //                   me, n_blocks, n_relevant_atoms*1e-3, distances_stored*1e-6, distances_stored*4e-6);
        if (echo > 2) std::printf("# rank#%i %d x %.3f k = %.6f M atom-block distances\n",
                              me, n_blocks, n_relevant_atoms*1e-3, distances_stored*1e-6);
        // if we knew the radii of compensation charges, core densities, valence densities, vbar, we could reduce this data item to a single bit, i.e. by 32x
        assert(n_relevant_atoms == atom_images.size());
        
        if (echo > 7) {
            std::printf("# rank#%i has contributing global_atom_ids {", me);
            for (auto const & ai : atom_images) {
                std::printf(" %i", ai.atom_id_);
            } // ai
            std::printf(" }, %ld in total\n", atom_images.size());
        } // echo

        return atom_images;
    } // get_neighborhood


    std::vector<uint32_t> find_unique_atoms(
          std::vector<atom_image_t> & atom_images
          , int const echo=0
    ) {
        // count how many periodic images each unique atom has
        std::map<int32_t,uint32_t> npi_map;
        for (auto const & ai : atom_images) { // 
            auto const global_atom_id = ai.atom_id_;
            ++npi_map[global_atom_id];
        } // ai
        // the size of this map determines how many different global_ids are involved
        auto const natoms = npi_map.size(); // number of contributing atoms
        if (echo > 0) std::printf("# found %ld different atoms in %ld atom images\n", natoms, atom_images.size());

        std::vector<uint32_t> global_atom_ids(natoms); // prepare function result
        std::map<int32_t,uint32_t> id_map; // this map translates from global_atom_ids to iatom-indices
        { // scope: 
            simple_stats::Stats<float> s(0);
            uint32_t iatom{0};
            for (auto const & atom : npi_map) { // loop over contributing atoms (not images)
                s.add(atom.second); // statistics about the number of periodic images
                auto const global_atom_id = atom.first;
                id_map[global_atom_id] = iatom;
                assert(global_atom_id >= 0);
                global_atom_ids[iatom] = global_atom_id;
                ++iatom;
            } // ia
            assert(natoms == iatom);
            auto const s_interval = s.interval();
            if (echo > 3) std::printf("# distribution of periodic images is %s\n", s_interval.c_str());
        } // scope

        // translate the members atom_id_ in atom_images from global ids into unique ids
        for (auto & ai : atom_images) {
            auto const global_atom_id = ai.atom_id_;
            auto const iatom = id_map[global_atom_id];
            ai.atom_id_ = iatom;
        } // ai

        return global_atom_ids;
    } // find_unique_atoms



    template <int PROJECT0_OR_ADD1, int n8=8>
    size_t sho_project0_or_add1(
          double      result[] // 0:grid_values[n8*n8*n8],  1:sho_coeff[nSHO(numax)], SHO-coefficients are zyx-ordered
        , double const input[] // 0:sho_coeff[nSHO(numax)], 1:grid_values[n8*n8*n8]
        , double const block_coords[3] // block origin
        , double const hg[3] // grid spacing
        , double const r_circum // == std::sqrt(pow2(h[0]) + pow2(h[1]) + pow2(h[2]))
        , double const atom_center[3] // where is the atom image
        , int const numax // SHO basis size
        , double const sigma // Gaussian spread
        , float const r_cut
        , int const echo=0 // log-level
    ) {
        if (1) { // check if the atom is too far from the center of the block
            auto const r2 = pow2(block_coords[0] + hg[0]*n8*.5 - atom_center[0])
                          + pow2(block_coords[1] + hg[1]*n8*.5 - atom_center[1])
                          + pow2(block_coords[2] + hg[2]*n8*.5 - atom_center[2]);
            if (r2 > pow2(r_cut + n8*.5*r_circum)) return 0;
        } // fast check

        float const r2cut = pow2(r_cut);
        assert(sigma > 0);
        auto const sigma_inv = 1./sigma;

        view3D<double> H1d(3,n8,sho_tools::n1HO(numax));
        float radius2[3][n8];
        for (int dir = 0; dir < 3; ++dir) {
            if (echo > 55) std::printf("\n# Hermite polynomials for %c-direction:\n", 'x' + dir);
            for (int i8 = 0; i8 < n8; ++i8) {
                auto const x = block_coords[dir] + (i8 + .5)*hg[dir] - atom_center[dir];
                Gauss_Hermite_polynomials(H1d(dir,i8), x*sigma_inv, numax);
                radius2[dir][i8] = pow2(x);
#ifdef    DEVEL
                if (echo > 55) {
                    std::printf("%g\t", x);
                    for (int nu = 0; nu <= numax; ++nu) {
                        std::printf(" %11.6f", H1d(dir,i8,nu));
                    } // nu
                    std::printf("\n");
                } // echo
#endif // DEVEL
            } // i8
        } // dir

        int const nSHO = sho_tools::nSHO(numax);

#ifdef    DEVEL
        if (1 == PROJECT0_OR_ADD1) {
            if (echo > 16) {
                std::printf("# addition coefficients ");
                for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                    std::printf(" %g", std::real(input[iSHO]));
                } // iSHO
                std::printf("\n");
            } // echo
        } // add
#endif // DEVEL

        size_t nhits{0};
        for (int iz{0}; iz < n8; ++iz) {
        for (int iy{0}; iy < n8; ++iy) {
        for (int ix{0}; ix < n8; ++ix) {

            if (radius2[0][ix] + radius2[1][iy] + radius2[2][iz] < r2cut) {

                int const izyx = (iz*n8 + iy)*n8 + ix;
                double val{0};
                if (0 == PROJECT0_OR_ADD1) {
                    val = input[izyx]; // load grid value
                } // project
                int iSHO{0};
                for (int nz = 0; nz <= numax; ++nz) {                    auto const H1d_z = H1d(2,iz,nz);
                    for (int ny = 0; ny <= numax - nz; ++ny) {           auto const H1d_y = H1d(1,iy,ny);
                        for (int nx = 0; nx <= numax - nz - ny; ++nx) {  auto const H1d_x = H1d(0,ix,nx);
                            auto const H3d = H1d_z * H1d_y * H1d_x;
                            if (1 == PROJECT0_OR_ADD1) {
                                val += input[iSHO] * H3d; // here, the addition happens
                            } else {
                                result[iSHO] += val * H3d; // here, the projection happens
                            }
                            ++iSHO; // in sho_tools::zyx_order
                        } // nx
                    } // ny
                } // nz
                assert( nSHO == iSHO );
                ++nhits;

//              if (echo > 6) std::printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                if (1 == PROJECT0_OR_ADD1) {
                    result[izyx] += val; // load-modify-store, must be atomic if threads are involved
                } // add

            } // inside radius

        }}} // ix iy iz

#ifdef    DEVEL
        if (0 == PROJECT0_OR_ADD1) {
            if (echo > 16) {
                std::printf("# projection coefficients ");
                for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                    std::printf(" %g", std::real(result[iSHO]));
                } // iSHO
                std::printf("\n");
            } // echo
        } // project
#endif // DEVEL

        return nhits;
    } // sho_project0_or_add1

    status_t add_to_grid(
          view2D<double> & grid_values // result: add to these grid array values, usually augmented_density
        , view2D<double> const & block_coords
        , uint32_t const n_blocks
        , data_list<double> const & atom_coeff // [natoms][nSHO]
        , std::vector<int32_t> const & lmax // lmax_qlm[natoms]
        , std::vector<double> const & sigma // sigma_cmp[natoms]
        , std::vector<atom_image_t> const & atom_images
        , double const hg[3] // == g.grid_spacings()
        , int const echo=0
    ) {
        auto const natoms = sigma.size();
        assert(natoms == lmax.size());
        assert(natoms == atom_coeff.nrows());
        auto const r_circum = std::sqrt(pow2(hg[0]) + pow2(hg[1]) + pow2(hg[2]));
        std::vector<size_t> hits_per_atom(sigma.size(), 0);
        for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // OMP PARALLEL
            size_t hits_per_block{0};
            for (auto const & ai : atom_images) {
                auto const iatom = ai.atom_id_;
                float const r_cut = 9*sigma[iatom];
                auto const hits = sho_project0_or_add1<1,8>(grid_values[ilb], atom_coeff[iatom], block_coords[ilb],
                                                            hg, r_circum, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                hits_per_atom[iatom] += hits;
                hits_per_block       += hits;
            } // ai
#ifdef    DEVEL
            if (echo > 17 && hits_per_block) std::printf("# %s: block #%i at %g %g %g Bohr has %.3f k hits\n", __func__,
                ilb, block_coords[ilb][0], block_coords[ilb][1], block_coords[ilb][2], hits_per_block*1e-3);
#endif // DEVEL
        } // ilb
#ifdef    DEVEL
        if (echo > 13) {
            for (size_t iatom{0}; iatom < sigma.size()*(echo > 0); ++iatom) {
                if (hits_per_atom[iatom]) std::printf("# %s: atom #%li has %.3f k hits\n", __func__, iatom, hits_per_atom[iatom]*1e-3);
            } // iatom
        } // echo
#endif // DEVEL
        return 0;
    } // add_to_grid

    status_t project_grid(
          data_list<double> & atom_coeff // result: atom_vlm[natoms][nSHO]
        , view2D<double> const & grid_values // project these grid array values, usually V_electrostatic
        , view2D<double> const & block_coords
        , uint32_t const n_blocks
        , std::vector<int32_t> const & lmax // lmax_vlm[natoms]
        , std::vector<double> const & sigma // sigma_cmp[natoms]
        , std::vector<atom_image_t> const & atom_images
        , double const hg[3] // == g.grid_spacings()
        , int const echo=0
    ) {
        auto const natoms = sigma.size();
        assert(natoms == lmax.size());
        assert(natoms == atom_coeff.nrows());
        auto const r_circum = std::sqrt(pow2(hg[0]) + pow2(hg[1]) + pow2(hg[2]));
        for (size_t iatom{0}; iatom < natoms; ++iatom) {
            auto const nSHO = sho_tools::nSHO(lmax[iatom]);
            set(atom_coeff[iatom], nSHO, 0.0); // clear
        } // iatom

        std::vector<size_t> hits_per_block(n_blocks, 0);
        for (auto const & ai : atom_images) { // OMP PARALLEL
            auto const iatom = ai.atom_id_;
            assert(iatom < natoms);
            float const r_cut = 9*sigma[iatom];
            size_t hits_per_atom{0};
            for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
                auto const hits = sho_project0_or_add1<0,8>(atom_coeff[iatom], grid_values[ilb], block_coords[ilb],
                                                            hg, r_circum, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                hits_per_block[ilb] += hits;
                hits_per_atom       += hits;
            } // ilb
#ifdef    DEVEL
            if (echo > 17 && hits_per_atom) std::printf("# %s: image of atom #%i at %g %g %g Bohr has %.3f k hits\n",
                                           __func__, iatom, ai.pos_[0], ai.pos_[1], ai.pos_[2], hits_per_atom*1e-3);
#endif // DEVEL
        } // ai
#ifdef    DEVEL
        if (echo > 13) {
            for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
                if (hits_per_block[ilb]) std::printf("# %s: block #%i at %g %g %g Bohr has %.3f k hits\n", __func__,
                    ilb, block_coords[ilb][0], block_coords[ilb][1], block_coords[ilb][2], hits_per_block[ilb]*1e-3);
            } // ilb
        } // echo
#endif // DEVEL
        return 0;
    } // project_grid


    template <typename real_t, int n8=8>
    double add_r2grid_to_block(
          real_t values[n8*n8*n8]
        , double const block_coords[3]
        , double const hg[3] // Cartesian grid spacings
        , double const r_circum // == std::sqrt(pow2(h[0]) + pow2(h[1]) + pow2(h[2]))
        , double const atom_center[3]
        , double const r2coeff[] // radial function on an r^2-grid
        , double const r_cut
        , double const factor=1
        , int const echo=0
        , int32_t const nr2=4096 // r^2-grid size
        , float const ar2=16.f // r^2-grid parameter
    ) {
        if (1) { // check if the atom is too far from the center of the block
            auto const r2 = pow2(block_coords[0] + hg[0]*n8*.5 - atom_center[0])
                          + pow2(block_coords[1] + hg[1]*n8*.5 - atom_center[1])
                          + pow2(block_coords[2] + hg[2]*n8*.5 - atom_center[2]);
            if (r2 > pow2(r_cut + n8*.5*r_circum)) return 0;
        } // fast check

        auto const r2cut = pow2(r_cut);
        double added_charge{0};
        size_t grid_points_inside{0};
        for (int iz{0}; iz < n8; ++iz) {
        for (int iy{0}; iy < n8; ++iy) {
        for (int ix{0}; ix < n8; ++ix) {
            auto const r2 = pow2(block_coords[0] + hg[0]*(ix + .5) - atom_center[0])
                          + pow2(block_coords[1] + hg[1]*(iy + .5) - atom_center[1])
                          + pow2(block_coords[2] + hg[2]*(iz + .5) - atom_center[2]);
            if (r2 < r2cut) {
                int const ir2 = int(ar2*r2);
                if (ir2 + 1 < nr2) {
                    double const w8 = ar2*r2 - ir2; // linear interpolation weight
                    auto const value_to_add = r2coeff[ir2] * (1 - w8) + r2coeff[ir2 + 1]*w8;
                    int const izyx = (iz*n8 + iy)*n8 + ix;
                    values[izyx] += factor*value_to_add;
                    added_charge += factor*value_to_add;
                }
                ++grid_points_inside;
            } // inside
        }}} // ix iy iz
#ifdef    DEVEL
        if (grid_points_inside && echo > 17) std::printf("# %ld grid points inside %g %s for grid block at [%g %g %g] Bohr\n",
            grid_points_inside, r_cut*Ang, _Ang, block_coords[0], block_coords[1], block_coords[2]);
#endif // DEVEL
        return added_charge;
    } // add_r2grid_to_block

    double add_r2grid_quantity(
          view2D<double> & grid_quantity // e.g. core_density
        , char const *const quantity_name
        , data_list<double> const & atom_r2coeff // e.g. atom_rhoc
        , std::vector<atom_image_t> const & atom_images
        , uint32_t const natoms
        , view2D<double> const & block_coords
        , uint32_t const n_blocks
        , real_space::grid_t const & g
        , MPI_Comm const comm
        , int const echo=0 // log level
        , double const factor=1
        , uint32_t const nr2=4096 // r^2-grid size
        , float const ar2=16.f // r^2-grid parameter
    ) {
        auto const me = mpi_parallel::rank(comm);
        // ToDo: get atom_r2coeff from atom owners
        auto const *const hg = g.grid_spacings();
        auto const r_circum = std::sqrt(pow2(hg[0]) + pow2(hg[1]) + pow2(hg[2]));

        float const r_cut = std::sqrt((nr2 - 1.)/ar2); // largest radius of the r^2-grid
        assert(ar2 > 0);
        assert(ar2 == 16.f);
        assert(ar2*pow2(r_cut) < nr2);
        assert(nr2 == 4096);

        double rho_total{0}; // prepare result
#ifdef    DEVEL
        std::vector<double> added_charge_per_atom(natoms, 0.0);
#endif // DEVEL
        for (auto const & ai : atom_images) {
            auto const iatom = ai.atom_id_;
            assert(0 <= iatom); assert(iatom < natoms);
#ifdef    DEVEL
            if (echo > 12) std::printf("# rank#%i adds to %s for atom#%i at position %g %g %g %s\n", 
                    me, quantity_name, iatom, ai.pos_[0]*Ang, ai.pos_[1]*Ang, ai.pos_[2]*Ang, _Ang);
#endif // DEVEL
            double added_charge{0};
            for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // local blocks
                added_charge += add_r2grid_to_block(grid_quantity[ilb], block_coords[ilb], hg, r_circum,
                                          ai.pos_, atom_r2coeff[iatom], r_cut, factor, echo, nr2, ar2);
            } // ilb
#ifdef    DEVEL
            if (echo > 13) std::printf("# rank#%i atom#%i %g electrons added for periodic image %i %i %i\n",
                                me, iatom, added_charge*g.dV(), ai.shifts_[0], ai.shifts_[1], ai.shifts_[2]);
            added_charge_per_atom[iatom] += added_charge;
#endif // DEVEL
            rho_total += added_charge;
        } // ia
        rho_total = mpi_parallel::sum(rho_total, comm)*g.dV();
#ifdef    DEVEL
        if (echo > 9) {
            double rho_total_checksum{0};
            for (int iatom{0}; iatom < natoms; ++iatom) {
                std::printf("# rank#%i atom#%i added %g electrons\n", me, iatom, added_charge_per_atom[iatom]*g.dV());
                rho_total_checksum += added_charge_per_atom[iatom];
            } // iatom
            std::printf("# rank#%i added %g electrons to %s\n", me, rho_total_checksum*g.dV(), quantity_name);
        } // echo
#endif // DEVEL
        if (echo > 3) std::printf("# added %g electrons to %s\n", rho_total, quantity_name);
        return rho_total;
    } // add_r2grid_quantity



    // =====================================================================    
    // ===== begin atom data communication routines ========================
    // =====================================================================    

    class AtomCommList_t {
    public:
        AtomCommList_t() {} // default constructor
        AtomCommList_t(
              size_t const n_all_atoms
            , std::vector<uint32_t> const & global_atom_ids // [natoms]
            , MPI_Comm const comm=MPI_COMM_WORLD
            , int const echo=0
        ) {
            comm_  = comm;
            auto const nprocs = mpi_parallel::size(comm); assert(nprocs > 0);
            auto const me     = mpi_parallel::rank(comm);
            int const na = (n_all_atoms + nprocs - 1 - me)/nprocs; // simple model, owner rank = global_atom_id % nprocs
            int const na_max =  (n_all_atoms + nprocs - 1)/nprocs;
            int const na_min =               (n_all_atoms)/nprocs;
            if (echo > 9) std::printf("# rank#%i has %d (min %d max %d) owned atoms\n", me, na, na_min, na_max);
            int const na_max8 = (na_max + 7) >> 3;
            if (echo > 5) std::printf("# %s: use MPI_Alltoall with %d--%d bits in %d Byte\n", __func__, na_min, na_max, na_max8);
            list_.resize(na);
            natoms_ = global_atom_ids.size();
#ifndef   HAS_NO_MPI
            view3D<uint8_t> bits(2, nprocs, na_max8, 0);
            auto bits_send = bits[0], bits_recv = bits[1];
            for (auto const & global_atom_id : global_atom_ids) {
                auto const atom_owner = global_atom_id % nprocs;
                auto const ia =         global_atom_id / nprocs;
             // bits(SEND,atom_owner,ia >> 3) |= (uint8_t(1) << (ia & 7));
                bits_send(atom_owner,ia >> 3) |= (uint8_t(1) << (ia & 7)); // set bit
            } // ai

            MPI_Alltoall(bits_send[0], na_max8, MPI_UINT8_T, bits_recv[0], na_max8, MPI_UINT8_T, comm);
            // Alternative: first use Alltoall for the number of atoms, then Alltoallv for the atom_ids or ias


            for (int ia{0}; ia < na; ++ia) { list_[ia].resize(0); } // init
            for (uint32_t rank{0}; rank < nprocs; ++rank) {
                for (int ia{0}; ia < na; ++ia) {
                    bool const atom_contributes = (bits_recv(rank,ia >> 3) >> (ia & 7)) & 1;
                    if (atom_contributes) {
                        list_[ia].push_back(rank);
                    }
                } // ia
                // all other bits must be unset
                for (int ia = na; ia < na_max8*8; ++ia) {
                    bool const atom_contributes = (bits_recv(rank,ia >> 3) >> (ia & 7)) & 1;
                    assert(!atom_contributes);
                } // ia
            } // rank
#else  // HAS_NO_MPI
            for (int ia{0}; ia < na; ++ia) { list_[ia].resize(1, 0); } // only rank zero contributes
#endif // HAS_NO_MPI
            
            if (echo > 9) {
                for (int ia{0}; ia < na; ++ia) {
                    auto const global_atom_id = ia*nprocs + me;
                    std::printf("# rank#%i communicates with %ld ranks for atom#%i global#%i\n", 
                                        me, list_[ia].size(), ia, global_atom_id);
                } // ia
            } // echo

        } // constructor

        std::vector<std::vector<uint32_t>> const & list() const { return list_; }
        MPI_Comm comm()   const { return comm_; }
        uint32_t natoms() const { return natoms_; }
    private:
        std::vector<std::vector<uint32_t>> list_;
        MPI_Comm comm_ = MPI_COMM_NULL;
        uint32_t natoms_ = 0;
    }; // class AtomCommList_t

    status_t atom_data_broadcast(
          data_list<double> & atom_data // result [natoms]
        , data_list<double> const & owner_data // input [na], only accessed in atom owner rank
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , int const echo=0 // log level
    ) {
        // send atom data updated by the atom owner to contributing MPI ranks
        status_t stat(0);
        auto const comm = atom_comm_list.comm();
        auto const nprocs = mpi_parallel::size(comm); assert(nprocs > 0);
        auto const me     = mpi_parallel::rank(comm);

        mpi_parallel::barrier(comm);

        // atom owners send their data
        auto const na = owner_data.nrows();
        auto const list = atom_comm_list.list();
        assert(na == list.size());
#ifndef   HAS_NO_MPI
        for (int ia{0}; ia < na; ++ia) { // loop over owned atoms
            auto const list_ia = list.at(ia);
            auto const count = owner_data.ncols(ia);
            for (auto const rank : list_ia) {
                if (rank != me) {
                    if (echo > 13) std::printf("# rank#%i %s: send %s, %d doubles for my owned atom#%i, global#%i to contributing rank#%i\n",
                                                       me, __func__, what, count, ia, ia*nprocs + me, rank);
                    MPI_Request send_request;
                    stat += MPI_Isend(owner_data[ia], count, MPI_DOUBLE, rank, ia, comm, &send_request);
                } // remote
            } // rank
        } // ia
        uint32_t irequest{0};
#else  // HAS_NO_MPI
        bool const remote_atom_is_error = (0 == control::get("mpi.fake.size", 0.));
#endif // HAS_NO_MPI


        // contributing atoms listen
        uint32_t const natoms = global_atom_ids.size();
        if (echo > 8) std::printf("# %s of %s, %d owned atoms to %d atoms\n", __func__, what, na, natoms);
        assert(natoms == atom_comm_list.natoms());
        std::vector<MPI_Request> recv_requests(natoms, MPI_REQUEST_NULL);
        for (uint32_t iatom{0}; iatom < natoms; ++iatom) { // loop over contributing atoms
            auto const global_atom_id = global_atom_ids[iatom];
            auto const atom_owner = global_atom_id % nprocs;
            auto const ia         = global_atom_id / nprocs;
            int const count = atom_data.ncols(iatom);
            if (atom_owner == me) {
                assert(count == owner_data.ncols(ia));
                if (echo > 11) std::printf("# rank#%i %s: copy %s, %d doubles for owned atom#%i to contributing atom#%i, global#%i, owner rank#%i\n",
                                                   me, __func__, what, count, ia, iatom, global_atom_id, atom_owner);
                set(atom_data[iatom], count, owner_data[ia]); // local copy
            } else {
#ifdef    HAS_NO_MPI
                if (remote_atom_is_error) error("cannot operate remote atoms without MPI, iatom= %i", iatom);
#else  // HAS_NO_MPI
                if (echo > 11) std::printf("# rank#%i %s: recv %s, %d doubles for owned atom#%i from owner rank#%i to contributing atom#%i, global#%i\n",
                                                   me, __func__, what, count, ia, atom_owner, iatom, global_atom_id);
                stat += MPI_Irecv(atom_data[iatom], count, MPI_DOUBLE, atom_owner, ia, comm, &recv_requests[irequest]);
                ++irequest;
#endif // HAS_NO_MPI
            }
        } // iatom

#ifndef   HAS_NO_MPI
        std::vector<MPI_Status> statuses(irequest);
        stat += MPI_Waitall(irequest, recv_requests.data(), statuses.data());
#endif // HAS_NO_MPI

        mpi_parallel::barrier(comm);

        if (stat) warn("failed for %s with status= %i", what, int(stat));
        return stat;
    } // atom_data_broadcast




    status_t atom_data_allreduce(
          data_list<double> & owner_data // result [na], only correct in atom owner rank
        , data_list<double> const & atom_data // input [natoms]
        , char const *const what
        , AtomCommList_t const & atom_comm_list
        , std::vector<uint32_t> const & global_atom_ids
        , double const factor=1 // scaling factor, usually g.dV(), the grid volume element
        , int const echo=0 // log level
    ) {
        // collect the unrenormalized atom_vlm data from contributing MPI ranks
        status_t stat(0);
        auto const comm = atom_comm_list.comm();
        auto const nprocs = mpi_parallel::size(comm);
        auto const me     = mpi_parallel::rank(comm);

        mpi_parallel::barrier(comm);

        // initialize the accumulators
        auto const na = owner_data.nrows();
        for (int ia{0}; ia < na; ++ia) {
            set(owner_data[ia], owner_data.ncols(ia), 0.0); // clear
        } // ia

        // contributing atoms
        uint32_t const natoms = global_atom_ids.size();
        if (echo > 8) std::printf("# %s of %s, %d atoms to %d owned atoms\n", __func__, what, natoms, na);
        assert(natoms == atom_comm_list.natoms());
        for (uint32_t iatom{0}; iatom < natoms; ++iatom) { // loop over contributing atoms
            auto const global_atom_id = global_atom_ids[iatom];
            auto const atom_owner = global_atom_id % nprocs;
            auto const ia         = global_atom_id / nprocs;
            int const count = atom_data.ncols(iatom);
            if (atom_owner == me) {
                assert(count == owner_data.ncols(ia));
                if (echo > 11) std::printf("# rank#%i %s:  add %s, %d doubles for owned atom#%i to contributing atom#%i, global#%i, owner rank#%i\n",
                                                   me, __func__, what, count, ia, iatom, global_atom_id, atom_owner);
                add_product(owner_data[ia], count, atom_data[iatom], factor); // local accumulation
            } else {
#ifdef    HAS_NO_MPI
                error("cannot operate remote atoms without MPI, iatom= %i", iatom);
#else  // HAS_NO_MPI
                if (echo > 11) std::printf("# rank#%i %s: send %s, %d doubles for contributing atom#%i, global#%i to owned atom#%i at owner rank#%i\n",
                                                   me, __func__, what, count, iatom, global_atom_id, ia, atom_owner);
                MPI_Request send_request;
                stat += MPI_Isend(atom_data[iatom], count, MPI_DOUBLE, atom_owner, ia, comm, &send_request);
#endif // HAS_NO_MPI
            }
        } // iatom

        // atom owners collect the data
        auto const list = atom_comm_list.list();
        assert(na == list.size());
#ifndef   HAS_NO_MPI
        std::vector<MPI_Request> recv_requests(1);
        for (int ia{0}; ia < na; ++ia) { // loop over owned atoms
            auto const list_ia = list.at(ia);
            auto const count = owner_data.ncols(ia);
            std::vector<double> contrib(count);
            for (auto const rank : list_ia) {
                if (rank != me) {
                    if (echo > 13) std::printf("# rank#%i %s: recv %s, %d doubles for my owned atom#%i, global#%i from contributing rank#%i\n",
                                                       me, __func__, what, count, ia, ia*nprocs + me, rank);
                    MPI_Status status;     // Mind that this is a blocking communication routine
                    stat += MPI_Recv(contrib.data(), count, MPI_DOUBLE, rank, ia, comm, &status);
                    add_product(owner_data[ia], count, contrib.data(), factor); // accumulation
                } // remote
            } // rank
        } // ia
#endif // HAS_NO_MPI

        mpi_parallel::barrier(comm);

        if (stat) warn("failed with status= %i", int(stat));
        return stat;
    } // atom_data_allreduce

    // =====================================================================    
    // ==== end of atom data communication routines ========================
    // =====================================================================    




    status_t SCF(int const echo) {
        status_t stat(0);
        SimpleTimer init_timer(strip_path(__FILE__), __LINE__, "init timer", 0);

        auto const comm = MPI_COMM_WORLD;
        auto const me = mpi_parallel::rank(comm);
        auto const nprocs = mpi_parallel::size(comm); assert(nprocs > 0);

        int const check = control::get("check", 0.); // check-mode

        bool needs_integrator{false};
        char const *const basis_method = control::get("basis", "Green-function"); // {Thomas-Fermi, Green-function}
        switch (*basis_method | 32) {
            case 't': break; // Thomas-Fermi model
            case 'g': needs_integrator = true; break; // Green-function model
            case 'n': warn("with +basis=%s --> none no new valence density is created!", basis_method); break;
            default:  error("not implemented +basis=%s", basis_method);
        } // switch

        // load geometry from +geometry.file=atoms.xyz
        real_space::grid_t g;    // entire grid descriptor
        view2D<double> xyzZ_all; // coordinates for all atoms
        int32_t n_all_atoms;     // number of all atoms (at most 2.147483647e9)

        { // scope: init_geometry_and_grid
            auto const stat_init = geometry_input::init_geometry_and_grid(g, xyzZ_all, n_all_atoms, 8, echo);
            if (stat_init) warn("init_geometry_and_grid returned status= %i", int(stat_init));
            stat += stat_init;
        } // scope

        // create a coarse grid descriptor
        real_space::grid_t gc(g[0] >> 1, g[1] >> 1, g[2] >> 1); // divide +grid.points by 2
        {
            assert(gc[0]*2 == g[0]); assert(gc[1]*2 == g[1]); assert(gc[2]*2 == g[2]); // g.grid_points must be an even number
            auto const gbc = g.boundary_conditions();
            gc.set_boundary_conditions(gbc);
            gc.set_cell_shape(g.cell, echo*0); // muted
            gc.set_grid_spacing(g.h[0]*2, g.h[1]*2, g.h[2]*2); // twice the grid spacing
            for (int d{0}; d < 3; ++d) { assert(0 == (gc[d] & 0x3)); } // all grid numbers must be a multiple of 4
            auto const max_grid_spacing = std::max(std::max(std::max(1e-9, gc.h[0]), gc.h[1]), gc.h[2]);
            if (echo > 1) std::printf("# use  %g %g %g  %s coarse grid spacing, corresponds to %.1f Ry\n",
                      gc.h[0]*Ang, gc.h[1]*Ang, gc.h[2]*Ang, _Ang, pow2(constants::pi/max_grid_spacing));
            g.set_boundary_conditions(boundary_condition::potential_bc(gbc[0]), 
                                      boundary_condition::potential_bc(gbc[1]), // map vacuum --> isolated, repeat --> periodic
                                      boundary_condition::potential_bc(gbc[2]));
        }
        double const grid_center[] = {g[0]*g.h[0]*.5, g[1]*g.h[1]*.5, g[2]*g.h[2]*.5}; // reference point for atomic positions

        // distribute the dense grid in 8x8x8 grid blocks to parallel owners
        parallel_poisson::load_balancing_t const lb(g, comm, 8, echo);
        if (echo > 1) { auto const nb = lb.grid_blocks(); std::printf("# use  %d %d %d  grid blocks\n", nb[0], nb[1], nb[2]); }

        parallel_poisson::parallel_grid_t const pg(g, lb, echo, "grid distribution");

        parallel_poisson::parallel_grid_t const pg_Interpolation(g, lb, echo, "Interpolation");


        // distribute the atom ownership:
        // simple distribution model: owner_rank == global_atom_id % nprocs, advantage: every rank can compute the owner
        //  
        // Example:  n_all_atoms=19, nprocs=7
        // rank      r#0   r#1   r#2   r#3   r#4   r#5   r#6
        //           a#0   a#1   a#2   a#3   a#4   a#5   a#6
        // owns      a#7   a#8   a#9   a#10  a#11  a#12  a#13
        //           a#14  a#15  a#16  a#17  a#18
        // na ==     3     3     3     3     3     2     2      na: number of atoms owned by this rank
        //
        // Alternative 1: instead of cyclic as above make something like owner_rank = (global_atom_id * nprocs)/n_all_atoms
        //
        // Alternative 2: make a weight grid pg.nb^3, init 0, add weight where atoms are, call load_balancer
        //     could be enriched by a workload model as function of Z
        //     would create some MPI closeness for systems with homogeneously distributed atoms

        if (echo > 0) std::printf("# distribute %.3f k atoms as %g atoms per rank on %d ranks\n", n_all_atoms*1e-3, n_all_atoms/double(nprocs), nprocs);
        int32_t const na = (n_all_atoms + nprocs - 1 - me)/nprocs; // simple model, owner_rank = global_atom_id % nprocs
        if (echo > 4) std::printf("# rank#%i has %d owned atoms\n", me, na);

        std::vector<double> Z_owned_atoms(na, 0.);
        for (int32_t ia{0}; ia < na; ++ia) {
            auto const gid = nprocs*ia + me;
            assert(0 <= gid); assert(gid < n_all_atoms);
            Z_owned_atoms[ia] = xyzZ_all(gid,3); // component 3 is the atomic number Z
        } // ia

        auto const n_blocks = pg.n_local();
        view2D<double> block_coords(n_blocks, 4, 0.0);
        auto atom_images_ = get_neighborhood(block_coords, n_all_atoms, xyzZ_all, pg, g, grid_center, echo);
        // for now atom_images_.atom_id_ are in [0, n_all_atoms)

        // ToDo: use get_neighborhood with r_cut + r_trunc to prefilter the relevant atomic images for the Green function method, so we don't have to pass xyzZ_all to them

        auto const global_atom_ids = find_unique_atoms(atom_images_, echo);
        auto const & atom_images = atom_images_;
        // from here atom_images.atom_id_ are in [0, global_atom_ids.size())



        uint32_t const natoms = global_atom_ids.size(); // number of locally contributing atoms
        if (echo > 4) std::printf("# rank#%i has %d locally contributing atoms\n", me, natoms);

        AtomCommList_t const atom_comm_list(n_all_atoms, global_atom_ids, comm, echo);

        float take_atomic_valence_densities{1}; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities
        if (echo > 2) std::printf("# take atomic valence densities with %g %%\n", take_atomic_valence_densities*100);

        char const *const pawdata_from = control::get("pawdata.from", "auto"); // 'a': auto generate, 'f': pawxml_import
        auto const pawdata_from_file = ('f' == (pawdata_from[0] | 32));
        if (echo > 2) std::printf("# use pawdata.from=%s  options {a, f} --> %s\n", pawdata_from, pawdata_from_file?"read from files":"generate");
        std::vector<int32_t> numax(na, pawdata_from_file ? -9 : -4); // -4: libliveatom, -9: load from pawxml files
        std::vector<int32_t> lmax_qlm(na, -1);      // expansion of qlm on owned atoms
        std::vector<int32_t> lmax_vlm(na, -1);      // expansion of vlm on owned atoms
        std::vector<int32_t> lmaxs_qlm(natoms, -1); // expansion of qlm on contributing atoms
        std::vector<int32_t> lmaxs_vlm(natoms, -1); // expansion of vlm on contributing atoms
        std::vector<int32_t> nr2(na,      1 << 12); // 4096 r^2-grid size on owned atoms
        std::vector<int32_t> nr2s(natoms, 1 << 12); // 4096 r^2-grid size on contributing atoms
        std::vector<float>   ar2(na,      16.f);    // r^2-grid parameter on owned atoms
        std::vector<float>   ar2s(natoms, 16.f);    // r^2-grid parameter on contributing atoms
        std::vector<double> sigma_cmp(na, 1.);      // spread of the Gaussian used in the compensation charges on owned atoms
        std::vector<double> sigmas_cmp(natoms, 1.); // spread of the Gaussian used in the compensation charges on contributing atoms

        // initialize and get sigma, lmax for each atom
        data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat; // for owned atoms
        data_list<double> atom_qzyx, atom_vzyx;                   // for owned atoms
        data_list<double> atoms_qzyx, atoms_vzyx;                 // for contributing atoms
        {
#ifdef    HAS_SINGLE_ATOM
            if (echo > 1) std::printf("# use single_atom::atom_update(what, na=%d, ...)\n", na);
#else  // HAS_SINGLE_ATOM
#ifdef    HAS_LIVE_ATOM
            // use linked library libliveatom
            if (echo > 1) std::printf("# use C-interface live_atom_update_(what, na=%d, ...)\n", na);
#else  // HAS_LIVE_ATOM
            if (echo > 1) std::printf("# missing live atom library -DHAS_LIVE_ATOM or -DHAS_SINGLE_ATOM\n");
#endif // HAS_LIVE_ATOM
#endif // HAS_SINGLE_ATOM

            std::vector<float> ionization(na, 0.f);
            stat += live_atom_update("initialize", na, Z_owned_atoms.data(), numax.data(), ionization.data(), (double**)1);
            stat += live_atom_update("lmax qlm",   na,    nullptr, lmax_qlm.data(), &take_atomic_valence_densities);
            stat += live_atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
            stat += live_atom_update("sigma cmp",  na, sigma_cmp.data());

            { // scope: initialized the data_lists for owned atoms
                view2D<uint32_t> num(6, na, 0); // how many
                for (int32_t ia{0}; ia < na; ++ia) {
                    num(0,ia) = pow2(1 + lmax_qlm.at(ia)); // qlm
                    num(1,ia) = pow2(1 + lmax_vlm.at(ia)); // vlm
                    auto const n_sho = sho_tools::nSHO(numax.at(ia)); assert(n_sho > 0);
                    num(2,ia) =   pow2(n_sho);             // aDm
                    num(3,ia) = 2*pow2(n_sho);             // aHm+aSm
                    num(4,ia) = sho_tools::nSHO(lmax_qlm[ia]); // number of coefficients to represent qlm compensation charges in a SHO basis
                    num(5,ia) = sho_tools::nSHO(lmax_vlm[ia]); // number of coefficients to represent vlm electrostatic projectors in a SHO basis
                } // ia
                // get memory in the form of data_list containers
                atom_qlm  = data_list<double>(na, num[0], 0.0); // charge multipole moments on owned atoms
                atom_vlm  = data_list<double>(na, num[1], 0.0); // electrostatic moments    on owned atoms
                atom_rho  = data_list<double>(na, num[2], 0.0); // atomic density matrices  on owned atoms
                atom_mat  = data_list<double>(na, num[3], 0.0); // atomic Hamiltonian       on owned atoms
                atom_qzyx = data_list<double>(na, num[4], 0.0); //                          on owned atoms
                atom_vzyx = data_list<double>(na, num[5], 0.0); //                          on owned atoms
            } // scope

            { // scope: group, broadcast and ungroup scalar atom data
                unsigned constexpr m8 = 8; // up to 8 scalars are transmitted as doubles
                std::vector<uint32_t> num(na, m8);
                data_list<double> atom_send(num, 0.0);
                for (int32_t ia{0}; ia < na; ++ia) {
                    atom_send(ia,0) = numax.at(ia);
                    atom_send(ia,1) = lmax_qlm.at(ia);
                    atom_send(ia,2) = lmax_vlm.at(ia);
                    atom_send(ia,3) = sigma_cmp.at(ia);
                    atom_send(ia,4) = nr2.at(ia);
                    atom_send(ia,5) = ar2.at(ia);
                    atom_send(ia,6) = ionization.at(ia);
                    atom_send(ia,7) = 0; // spare
                } // ia

                num.resize(natoms, m8);
                data_list<double> atoms_recv(num, 0.0);

                stat += atom_data_broadcast(atoms_recv, atom_send, "eight atom scalars", atom_comm_list, global_atom_ids, echo);

                for (uint32_t iatom{0}; iatom < natoms; ++iatom) {
                    lmaxs_qlm.at(iatom)  = atoms_recv(iatom,1);
                    lmaxs_vlm.at(iatom)  = atoms_recv(iatom,2);
                    sigmas_cmp.at(iatom) = atoms_recv(iatom,3); 
                    nr2s.at(iatom)       = atoms_recv(iatom,4);
                    ar2s.at(iatom)       = atoms_recv(iatom,5);
                } // iatom
            } // scope

            { // scope: initialized the data_lists for addition and projection with contributing atoms
                view2D<uint32_t> num(2, natoms, 0); // how many
                for (uint32_t iatom{0}; iatom < natoms; ++iatom) {
                    num(0,iatom) = sho_tools::nSHO(lmaxs_qlm[iatom]);
                    num(1,iatom) = sho_tools::nSHO(lmaxs_vlm[iatom]);
                } // iatom
                atoms_qzyx = data_list<double>(natoms, num[0], 0.0); // charge coefficients  on contributing atoms
                atoms_vzyx = data_list<double>(natoms, num[1], 0.0); // electrostatic coeffs on contributing atoms
            } // scope

        } // scope

        // the r^2-grids are used to bring radial quantities to a Cartesian grid
        data_list<double> atom_vbar(nr2, 0.0); // zero potentials on owned atoms
        data_list<double> atom_rhoc(nr2, 0.0); // core densities  on owned atoms

        data_list<double> atoms_vbar(nr2s, 0.0); // zero potentials on contributing atoms
        data_list<double> atoms_rhoc(nr2s, 0.0); // core densities  on contributing atoms



        std::vector<int32_t> numax_prj;
        std::vector<double>  sigma_prj;
        energy_contour::Integrator integrator;
        if (needs_integrator) {
            if (echo > 0) std::printf("\n# Initialize energy contour integrator\n");

            // Mind: currently the Green function method requires all atoms
            //     as it has to tell apart atomic images from atomic copies
            std::vector<double> xyzZinso(0);
            { // scope: determine additional info
                numax_prj.resize(na, 0);
                sigma_prj.resize(na, 1.);
                stat += live_atom_update("projectors", na, sigma_prj.data(), numax_prj.data());

                view2D<double> numax_sigma(n_all_atoms, 2, 0.0);
                for (int32_t ia{0}; ia < na; ++ia) { // loop over owned atoms
                    assert(numax_prj.at(ia) == numax.at(ia) && "inconsist between 'projectors' and 'initialize' call");
                    auto const gid = nprocs*ia + me; // global_atom_id
                    assert(0 <= gid); assert(gid < n_all_atoms);
                    numax_sigma(gid,0) = numax_prj.at(ia);
                    numax_sigma(gid,1) = sigma_prj.at(ia);
                } // ia
                mpi_parallel::max(numax_sigma[0], n_all_atoms*2, comm); // this is potentially slow

                xyzZinso.resize(8*n_all_atoms);
                for (int32_t gid{0}; gid < n_all_atoms; ++gid) { // another loop over all atoms, TODO can we avoid this?
                    set(&xyzZinso[gid*8], 4, xyzZ_all[gid]); // copy position and atomic number
                    xyzZinso[gid*8 + 4] = gid;
                    xyzZinso[gid*8 + 5] = numax_sigma(gid,0);
                    xyzZinso[gid*8 + 6] = numax_sigma(gid,1);
                    xyzZinso[gid*8 + 7] = 0; // spare
                } // gid
            } // scope

            // envoke the constructor energy_contour::Integrator
            integrator = energy_contour::Integrator(gc, xyzZinso, echo, check);

            // setup communication infrastructure for atom_mat
            auto const & target_global_atom_ids = integrator.plan_->dyadic_plan.global_atom_ids;
            std::vector<int64_t> owned_global_atom_ids(na);
            for (int ia{0}; ia < na; ++ia) {
                owned_global_atom_ids[ia] = nprocs*ia + me;
            } // ia
            assert(0 == (xyzZinso.size() & 0x7)); // make sure it is divisible by 8
            uint32_t const n_all_atoms = xyzZinso.size() >> 3; // divide by 8
            std::vector<green_parallel::rank_int_t> atom_owner_rank(n_all_atoms, 0);
            auto const nprocs = mpi_parallel::size(); // comm=MPI_COMM_WORLD?
            for (uint32_t gid{0}; gid < n_all_atoms; ++gid) {
                atom_owner_rank[gid] = gid % nprocs;
            } // gid
            uint32_t const nb[] = {n_all_atoms, 0, 0};
            integrator.plan_->matrices_requests = green_parallel::RequestList_t(
                target_global_atom_ids, owned_global_atom_ids, atom_owner_rank.data(), nb, echo, "atom matrices");
            if (echo > 1) std::printf("\n");
        } // needs_integrator

        xyzZ_all = view2D<double>(0, 0, 0.0); // clear, xyzZ_all should not be used after this




        // allocate CPU memory for grid arrays
        view2D<double> augmented_density(n_blocks, 8*8*8, 0.0);
        view2D<double> V_electrostatic(  n_blocks, 8*8*8, 0.0);
        view2D<double> core_density(     n_blocks, 8*8*8, 0.0);
        view3D<double> valence_density(2,n_blocks, 8*8*8, 0.0); // 0: density, 1: energy derivative w.r.t. the Fermi level


        double constexpr Y00 = .28209479177387817; // == 1/sqrt(4*pi)
        double constexpr Y00sq = pow2(Y00);        // == 1/(4*pi)

        // configure the electrostatic solver == Poisson solver
        auto  const es_method    = control::get("poisson.method", "cg"); // {cg, sd}
        int   const es_echo      = control::get("poisson.echo", echo/2);
        float const es_threshold = control::get("poisson.threshold", 3e-8);
        int   const es_maxiter   = control::get("poisson.maxiter", 200.);
        int   const es_miniter   = control::get("poisson.miniter", 3.);
        int   const es_restart   = control::get("poisson.restart", 4096.);

        // configure the self-consistency loop
        int   const scf_maxiter  = control::get("scf.maxiter", 1.);
        if (echo > 2) std::printf("#\n# Start up to %d SCF-iterations\n#\n", scf_maxiter);

        double E_Fermi{0}; // Fermi level

        double nve{0}; // non-const total number of valence electrons
        { // scope: determine the number of valence electrons
            auto const keyword_valence_electrons = "valence.electrons";
            if ('a' == (*control::get(keyword_valence_electrons, "auto") | 32)) { // | 32 converts to lowercase
                std::vector<double> n_electrons_a(na, 0.); // number of valence electrons added by each owned atom
                stat += live_atom_update("#valence electrons", na, n_electrons_a.data());
                nve = std::accumulate(n_electrons_a.begin(), n_electrons_a.end(), 0.0);
                nve = mpi_parallel::sum(nve, comm);
                if (echo > 0) std::printf("\n# %s=auto --> %g valence electrons\n\n", keyword_valence_electrons, nve);
                control::set(keyword_valence_electrons, nve, control::echo_set_without_warning);
            } else {
                nve = control::get(keyword_valence_electrons, 0.0);
                if (echo > 0) std::printf("\n# %s=%g\n\n", keyword_valence_electrons, nve);
            } // auto
        } // scope
        double const n_valence_electrons = nve; // do not use nve beyond this point

        sho_unitary::Unitary_SHO_Transform const unitary(9);

        {
            auto const init_time = init_timer.stop();
            if (echo > 2) { std::printf("# SCF initialization took %g seconds\n", init_time); }
        }

        // energy contributions
        double grid_kinetic_energy{0}, grid_electrostatic_energy{0}, grid_xc_energy{0}, total_energy{0};

        simple_stats::Stats<> scf_iteration_times, poisson_solver_times, green_function_times;

        int scf_iteration{0};
        bool scf_run{true};
        do { // self-consistency loop

            ++scf_iteration;
            if (echo > 3) std::printf("#\n# Start SCF-iteration#%i\n#\n", scf_iteration);
            mpi_parallel::barrier(comm);

            char scf_iteration_label[64];
            std::snprintf(scf_iteration_label, 64, "SCF-iteration#%i", scf_iteration);
            SimpleTimer scf_iteration_timer(strip_path(__FILE__), __LINE__, scf_iteration_label, echo*(0 == check));


            // get smooth core densities
            stat += live_atom_update("core densities", na, 0, nr2.data(), 0, atom_rhoc.data());

            if (take_atomic_valence_densities > 0) {
                auto & atom_rhov = atom_vbar; // use memory of vbar for the moment
                // get smooth spherical valence density
                stat += live_atom_update("valence densities", na, 0, nr2.data(), 0, atom_rhov.data());
                for (int32_t ia{0}; ia < na; ++ia) { // loop over owned atoms
                    // add valence density to r^2-gridded core density
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core electrons\n",         me, ia, integrate_r2grid(atom_rhoc[ia], nr2[ia]));
                    add_product(atom_rhoc[ia], nr2[ia], atom_rhov[ia], take_atomic_valence_densities*1.);
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core+valence electrons\n", me, ia, integrate_r2grid(atom_rhoc[ia], nr2[ia]));
                } // ia
            } // take_atomic_valence_densities > 0
            stat += atom_data_broadcast(atoms_rhoc, atom_rhoc, "core densities", atom_comm_list, global_atom_ids, echo);

            auto const total_charge_added = add_r2grid_quantity(core_density, "smooth core density", atoms_rhoc,
                                            atom_images, natoms, block_coords, n_blocks, g, comm, echo, Y00sq);
            if (echo > 0) std::printf("# %g electrons added as smooth core density\n", total_charge_added);

            // compose density
            set(augmented_density[0], n_blocks*size_t(8*8*8), core_density[0]);
            if (take_atomic_valence_densities < 1) { // less than 100% atomic valence densities
                add_product(augmented_density[0], n_blocks*size_t(8*8*8), valence_density(0,0), 1.);
            } // take_atomic_valence_densities < 1

            print_stats(augmented_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# smooth density");

            view2D<double> V_xc(n_blocks, 8*8*8, 0.0);
            { // scope: eval the XC potential and energy on the dense grid
                double E_xc{0}, E_dc{0};
                auto const *const density = augmented_density[0];
                auto       *const potential = V_xc[0];
                // double rho_max{0}; int64_t i_max{-1};
                for (size_t i = 0; i < n_blocks*size_t(8*8*8); ++i) {
                    auto const rho_i = density[i];
                    // if (rho_i > rho_max) { rho_max = rho_i; i_max = i; }
                    double vxc_i;
                    auto const exc_i = exchange_correlation::LDA_kernel(rho_i, vxc_i);
                    E_xc += rho_i*exc_i;
                    E_dc += rho_i*vxc_i; // double counting correction
                    potential[i] = vxc_i;
                } // i
             // if (echo > 2) std::printf("# rank#%i rho_max= %g a.u. at index= %lli\n", me, rho_max, i_max);
                E_xc = mpi_parallel::sum(E_xc, comm) * g.dV(); // scale with volume element
                E_dc = mpi_parallel::sum(E_dc, comm) * g.dV(); // scale with volume element
                if (echo > 2) std::printf("# exchange-correlation energy on grid %.9f %s, double counting %.9f %s\n", E_xc*eV, _eV, E_dc*eV, _eV);
                grid_xc_energy = E_xc;
            } // scope
            print_stats(V_xc[0], n_blocks*size_t(8*8*8), comm, echo > 1, 0, "# smooth exchange-correlation potential", eV, _eV);

            stat += live_atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());

            for (int ia{0}; ia < na; ++ia) {
                auto const global_atom_id = ia*nprocs + me;
                auto const stat_den = sho_projection::denormalize_electrostatics(atom_qzyx[ia], atom_qlm[ia], lmax_qlm[ia], sigma_cmp[ia], unitary, echo);
                if (stat_den) warn("denormalize_electrostatics failed with status= %i for atom#%i", int(stat_den), global_atom_id);
                stat += stat_den;
            } // ia
            stat += atom_data_broadcast(atoms_qzyx, atom_qzyx, "compensator multipole moments", atom_comm_list, global_atom_ids, echo);

            add_to_grid(augmented_density, block_coords, n_blocks, atoms_qzyx, lmaxs_qlm, sigmas_cmp, atom_images, g.grid_spacings(), echo);

            print_stats(augmented_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# smooth augmented_density");


            // ====================================================================================
            // ====================================================================================
            // =====                    Poisson equation                    =======================
            // ====================================================================================
            if (0 == check) { // scope: Poisson equation
                if (echo > 3) std::printf("#\n# Solve the Poisson equation iteratively with %d ranks in SCF-iteration#%i\n#\n", nprocs, scf_iteration);
                std::snprintf(scf_iteration_label, 64, "Poisson in SCF-iteration#%i", scf_iteration);
                SimpleTimer poisson_timer(strip_path(__FILE__), __LINE__, scf_iteration_label, echo);
                double E_es{0};
                float es_residual{0};

                auto const es_stat = parallel_poisson::solve(V_electrostatic[0], augmented_density[0],
                                                pg, *es_method, es_echo, es_threshold,
                                                &es_residual, es_maxiter, es_miniter, es_restart, &E_es);

                poisson_solver_times.add(poisson_timer.stop());
                if (echo > 2) std::printf("# Poisson equation %s, residual= %.2e a.u.\n#\n", es_stat?"failed":"converged", es_residual);
                stat += es_stat;
                if (es_stat && 0 == me) warn("parallel_poisson::solve returned status= %i", int(es_stat));
                grid_electrostatic_energy = 0.5*E_es; // store
                if (echo > 3) std::printf("# smooth electrostatic grid energy %.9f %s\n", grid_electrostatic_energy*eV, _eV);
            } else {
                if (echo > 0) std::printf("\n# skip Poisson equation for the electrostatic potential due to +check=%d\n\n", check);
                set(V_electrostatic[0], n_blocks*size_t(8*8*8), 0.0);
            }
            // ====================================================================================

            print_stats(V_electrostatic[0], n_blocks*size_t(8*8*8), comm, echo > 0, 0, "# smooth electrostatic potential", eV, _eV);

            // project the electrostatic grid onto the localized compensation charges
            project_grid(atoms_vzyx, V_electrostatic, block_coords, n_blocks, lmaxs_vlm, sigmas_cmp, atom_images, g.grid_spacings(), echo);

            stat += atom_data_allreduce(atom_vzyx, atoms_vzyx, "projected electrostatic potential", atom_comm_list, global_atom_ids, g.dV(), echo);
            for (int32_t ia{0}; ia < na; ++ia) {
                auto const global_atom_id = ia*nprocs + me;
                auto const stat_ren = sho_projection::renormalize_electrostatics(atom_vlm[ia], atom_vzyx[ia], lmax_vlm[ia], sigma_cmp[ia], unitary, echo);
                if (stat_ren) warn("renormalize_electrostatics failed with status= %i for atom#%i", int(stat_ren), global_atom_id);
                stat += stat_ren;
#ifdef    DEVEL
                if (echo > 7) {
                    std::printf("# potential projection for atom #%i v_00 = %.9f %s\n", global_atom_id, atom_vlm[ia][00]*Y00*eV, _eV);
                    int const ellmax_show = std::min(lmax_vlm[ia], 2);
                    for (int ell = 1; ell <= ellmax_show; ++ell) {
                        double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                        int const ilm0 = sho_tools::lm_index(ell, -ell);
                        std::printf("# potential projection for atom #%i v_%im =", global_atom_id, ell);
                        printf_vector(" %.6f", &atom_vlm[ia][ilm0], 2*ell + 1, nullptr, unitfactor);
                        std::printf(" %s %s^%i\n", _eV, _Ang, -ell);
                    } // ell
                } // echo
#endif // DEVEL
            } // ia

            // compose and coarsen total effective potential

            view2D<double> V_effective(n_blocks, 8*8*8, 0.0);
            { // scope: add potentials
                set(        V_effective[0], n_blocks*size_t(512), V_electrostatic[0]);
                add_product(V_effective[0], n_blocks*size_t(512), V_xc[0], 1.);
            } // scope
            print_stats(V_effective[0], n_blocks*size_t(8*8*8), comm, echo > 0, 0, "# smooth effective potential", eV, _eV);

            float potential_mixing_ratio[] = {.5}; // {potential}
            stat += live_atom_update("update", na, 0, 0, potential_mixing_ratio, atom_vlm.data());
            stat += live_atom_update("hamiltonian", na, 0, 0, 0, atom_mat.data());
            stat += live_atom_update("zero potentials", na, 0, nr2.data(), 0, atom_vbar.data());

            stat += atom_data_broadcast(atoms_vbar, atom_vbar, "zero potentials", atom_comm_list, global_atom_ids, echo); // assume vbar is updated every scf-iteration

            add_r2grid_quantity(V_effective, "smooth effective potential", atoms_vbar,
                                atom_images, natoms, block_coords, n_blocks, g, comm, echo*0, Y00);

            print_stats(V_effective[0], n_blocks*size_t(8*8*8), comm, echo > 0, 0, "# smooth effective potential", eV, _eV);


            view2D<double> new_valence_density(n_blocks, 8*8*8, 0.0);

            switch (*basis_method | 32) {
            case 't':
            {
                if (echo > 0) std::printf("# +basis=%s --> Thomas-Fermi model\n", basis_method);
                // apply Thomas-Fermi approximation
                auto const stat_TF = new_density_Thomas_Fermi(new_valence_density[0], E_Fermi, V_effective[0],
                                n_blocks*size_t(512), comm, n_valence_electrons, g.dV(), echo);
                stat += stat_TF;
                if (stat_TF && 0 == me) warn("# new_density_Thomas_Fermi returned status= %i", int(stat_TF));
                print_stats(new_valence_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# new Thomas-Fermi density");
            }
            break;

            case 'g':
            {
                std::snprintf(scf_iteration_label, 64, "Green function in SCF-iteration#%i", scf_iteration);
                SimpleTimer green_timer(strip_path(__FILE__), __LINE__, scf_iteration_label, 0);

                if (echo > 0) std::printf("# +basis=%s --> Green-function model\n", basis_method);
                view2D<double> V_coarse(n_blocks, 4*4*4, 0.0);
                for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // parallel loop over local blocks
                    block_average(V_coarse[ilb], V_effective[ilb]);
                } // ilb
                print_stats(V_coarse[0], n_blocks*size_t(4*4*4), comm, echo > 0, 0, "# coarse effective potential", eV, _eV);

                double band_bottom{-1.};
                { // scope: extract the highest core state energy and an estimate for the lowest valence state energy
                    view2D<float> extreme_energy_a(3, na, 0.); // extremal energy for {core, semicore, valence}
                    stat += live_atom_update("#core electrons"   , na, 0, 0, extreme_energy_a[0]);
                    // ToDo: semicore states ...
                    stat += live_atom_update("#valence electrons", na, 0, 0, extreme_energy_a[2]);
                    double extreme[] = {-9e9, 0.0, 9e9};
                    for (int32_t ia{0}; ia < na; ++ia) {
                        extreme[0] = std::max(extreme[0], 1.*extreme_energy_a(0,ia));
                     // extreme[1] =          extreme[1]  +  extreme_energy_a(1,ia) ; 
                        extreme[2] = std::min(extreme[2], 1.*extreme_energy_a(2,ia));
                    } // ia
#ifdef    DEBUGGPU
                    if (echo > 5) std::printf("# rank#%i has energy extrema at %g and %g %s\n", me, extreme[0]*eV, extreme[2]*eV, _eV);
#endif // DEBUGGPU
                    extreme[0] = mpi_parallel::max(extreme[0], comm); // global maximum
                 // extreme[1] = mpi_parallel::sum(extreme[1], comm)/std::max(n_all_atoms, 1);
                    extreme[2] = mpi_parallel::min(extreme[2], comm); // global minimum
                    // ToDo: semicore contour ...
                    if (extreme[2] < 8e9) {
                        if (echo > 2) std::printf("# energy contour estimates lowest valence state at %g %s\n", extreme[2]*eV, _eV);
                        if (0.0 == E_Fermi) { // Fermi level is exactly 0.0 when it is initialized
                            E_Fermi = extreme[2];
                            if (echo > 0) std::printf("# initialize Fermi level as %g %s\n", E_Fermi*eV, _eV);
                        }
                        if (extreme[0] > -8e9) {
                            if (echo > 2) std::printf("# energy contour assumes highest core state at %g %s\n", extreme[0]*eV, _eV);
                            auto const suggest_bottom = 0.25*extreme[0] + 0.75*extreme[2];
                            band_bottom = suggest_bottom - E_Fermi;
                            if (echo > 0) std::printf("# suggest band.bottom %g %s below the Fermi level\n", band_bottom*eV, _eV);
                        }
                    }
                    control::set("energy_contour.band.bottom", band_bottom); // avoid changing the interface for now
                } // scope

                // call energy-contour integration to find a new density
                auto const stat_Gf = integrator.integrate(new_valence_density[0], E_Fermi, V_coarse[0], atom_mat, numax_prj, sigma_prj,
                                                          pg_Interpolation, n_valence_electrons, g.dV(), echo, check);
                stat += stat_Gf;
                if (stat_Gf && 0 == me) warn("# energy_contour::integration returned status= %i", int(stat_Gf));

                auto const green_function_took = green_timer.stop();
                mpi_parallel::barrier(comm); // wait until other ranks have finished the Green function solution
                simple_stats::Stats<> green_time_stats;
                green_time_stats.add(green_function_took);
                mpi_parallel::allreduce(green_time_stats);
                if (0 == check && echo > 2) {
                    std::printf("# Green function solution in SCF-iteration#%i took %s seconds\n", 
                        scf_iteration, green_time_stats.interval().c_str());
                } // echo
                green_function_times.add(green_time_stats.max());
            }
            break;

            case 'n':
            {
                if (echo > 0) std::printf("# +basis=%s --> none, valence density is not updated!\n", basis_method);
            }
            break;

            default: error("not implemented +basis=%s", basis_method);
            } // switch basis_method

            print_stats(new_valence_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# new valence density");


            auto const E_dcc = dot_product(n_blocks*size_t(8*8*8), new_valence_density[0], V_effective[0]);
            double const double_counting_correction = mpi_parallel::sum(E_dcc, comm) * g.dV();
            if (echo > 1) std::printf("\n# grid double counting %.9f %s\n\n", double_counting_correction*eV, _eV);

            grid_kinetic_energy = E_Fermi - double_counting_correction;
            if (echo > 1) std::printf("\n# grid kinetic energy %.9f %s (take %.1f %%)\n\n",
                    grid_kinetic_energy*eV, _eV, (1. - take_atomic_valence_densities)*100);

            float rho_mixing_ratios[] = {.5, .5, .5}; // for spherical {core, semicore, valence} density
            stat += live_atom_update("atomic density matrices", na, 0, 0, rho_mixing_ratios, atom_rho.data());

            double const mix_new = rho_mixing_ratios[2], mix_old = 1. - mix_new;
            scale(valence_density(0,0),       n_blocks*size_t(8*8*8), mix_old);
            add_product(valence_density(0,0), n_blocks*size_t(8*8*8), new_valence_density[0], mix_new);


            // compute the total energy
            std::vector<double> atomic_energy_diff(na, 0.0);

            bool const total_energy_details = true;
            if (total_energy_details) {
                int const nE = align<1>(energy_contribution::max_number);
                std::vector<int32_t> nEa(na, nE);
                data_list<double> atom_contrib(nEa, 0.0);
                stat += live_atom_update("energies", na, atomic_energy_diff.data(), 0, 0, atom_contrib.data());
                std::vector<double> Ea(nE, 0.0);
                for (int32_t ia = 0; ia < na; ++ia) {
                    add_product(Ea.data(), nE, atom_contrib[ia], 1.); // all atomic weight factors are 1.0
                } // ia
                mpi_parallel::sum(Ea.data(), nE, comm);
                if (echo > 7) std::printf("\n# sum of atomic energy contributions without grid contributions:\n");
                energy_contribution::show(Ea.data(), echo - 7, eV, _eV);
                // now add grid contributions
                Ea[energy_contribution::KINETIC] += grid_kinetic_energy * (1. - take_atomic_valence_densities);
                Ea[energy_contribution::EXCHANGE_CORRELATION] += grid_xc_energy;
                Ea[energy_contribution::ELECTROSTATIC] += grid_electrostatic_energy;
                // reconstruct total energy from its contributions KINETIC + ES + XC
                Ea[energy_contribution::TOTAL] = Ea[energy_contribution::KINETIC]
                                               + Ea[energy_contribution::ELECTROSTATIC]
                                               + Ea[energy_contribution::EXCHANGE_CORRELATION];
                if (echo > 7) std::printf("\n# energy contributions including grid contributions:\n");
                energy_contribution::show(Ea.data(), echo - 3, eV, _eV);
            } else {
                stat += live_atom_update("energies", na, atomic_energy_diff.data());
            } // total_energy_details

            double atomic_energy_corrections{0};
            for (int32_t ia = 0; ia < na; ++ia) {
                atomic_energy_corrections += atomic_energy_diff[ia];
            } // ia
            atomic_energy_corrections = mpi_parallel::sum(atomic_energy_corrections, comm);

            total_energy = grid_kinetic_energy * (1. - take_atomic_valence_densities)
                         + grid_xc_energy
                         + grid_electrostatic_energy
                         + atomic_energy_corrections;
            if (echo > 0) { std::printf("\n# total energy %.9f %s\n\n", total_energy*eV, _eV); std::fflush(stdout); }



            scf_run = (scf_iteration < scf_maxiter && 0 == check); // run only 1 iteration in check mode (check==1)

            mpi_parallel::barrier(comm);

            scf_iteration_times.add(scf_iteration_timer.stop());

        } while (scf_run); // self-consistency loop 


        if (echo > 1) { std::printf("\n# SCF-iterations stopped after %d of %d iterations\n", scf_iteration, scf_maxiter); }
        if (0 == check) {
            if (echo > 3) {
                auto const & s = poisson_solver_times;
                std::printf("# %ld Poisson solutions in %s took %g seconds\n", s.tim(), s.interval().c_str(), s.sum());
            } // echo
            if (echo > 2 && green_function_times.tim() > 0) {
                auto const & s = green_function_times;
                std::printf("# %ld Green functions in %s took %g seconds\n", s.tim(), s.interval().c_str(), s.sum());
            } // echo
            if (echo > 2) {
                auto const & s = scf_iteration_times;
                std::printf("# %ld SCF-iterations in %s took %g seconds\n", s.tim(), s.interval().c_str(), s.sum());
            } // echo
        } // do not report execution times of the check mode


        if (1) { // scope: report GPU memory consumption
            int const memory_show = control::get("green_memory.show", 0.);
            if (memory_show) std::printf("# rank#%i green_memory now= %.9f GByte, max= %.9f GByte\n",
                    me, green_memory::total_memory_now()*1e-9, green_memory::high_water_mark()*1e-9);
        } // scope

        stat += live_atom_update("memory cleanup", na);
        return stat;
    } // SCF



















#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

    status_t test_r2grid_integrator(int const echo=0, int const nr2=4096, float const ar2=16) {
        status_t stat(0);
        std::vector<double> vec(nr2, 1.0);
        { // integrate 1.0 with r^2*dr up to the outer max radius R
            auto const R = std::sqrt((nr2 - 1.)/ar2);
            auto const vol = integrate_r2grid(vec.data(), nr2, ar2), ref = pow3(R)/3, dev = vol - ref;
            if (echo > 3) std::printf("# %s: found %g expect %g dev= %.3f %%\n", __func__, vol, ref, dev/ref*100);
            stat += (std::abs(dev) > ref*2e-4);
        }
        { // prepare a Gaussian
            for (int ir2{0}; ir2 < nr2; ++ir2) { vec[ir2] = std::exp(-ir2*(.125/ar2)); }
            auto const vol = integrate_r2grid(vec.data(), nr2, ar2), ref = std::sqrt(constants::pi*32), dev = vol - ref;
            if (echo > 3) std::printf("# %s: found %g expect %g dev= %.3f %%\n", __func__, vol, ref, dev/ref*100);
            stat += (std::abs(dev) > ref*2e-4);
        }
        return stat;
    } // test_r2grid_integrator

    status_t test_scf(int const echo=0) {
        status_t stat(0);
        stat += SCF(echo);
        return stat;
    } // test_scf

    status_t all_tests(int const echo) {
        status_t stat(0);
        stat += test_r2grid_integrator(echo);
        auto const already_initialized = mpi_parallel::init();
        int64_t const echo_mask = control::get("verbosity.mpi.mask", 1.); // -1:all, 0:no_one, 1:master only, 5:rank#0 and rank#2, ...
        auto const myrank = mpi_parallel::rank();
        int const verbose_rank = (-1 == echo_mask) ? 1 : ((myrank < 53)*((echo_mask >> myrank) & 1));
        auto const echo_rank = verbose_rank*mpi_parallel::max(echo);
        stat += test_scf(echo_rank);
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_potential
