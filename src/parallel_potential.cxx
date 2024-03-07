// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <cstdint> // int32_t, uint32_t
#include <cmath> // std::exp, ::sqrt, ::abs
#include <map> // std::map<Key,T>

#include "parallel_potential.hxx"

// #include "single_atom.h" // live_atom_##NAME##_ // use the libliveatom library
#include "single_atom.hxx" // ::atom_update

#include "status.hxx" // status_t
#include "control.hxx" // ::get
#include "real_space.hxx" // ::grid_t
#include "display_units.h" // Ang, _Ang
#include "self_consistency.hxx" // ::init_geometry_and_grid
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

namespace parallel_potential {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

    inline status_t live_atom_update(
        char const *const what    // selector string
      , int32_t  const natoms     // number of atoms
      , double  *const dp=nullptr // quantities (input/output) double  dp[natoms]
      , int32_t *const ip=nullptr // quantities (input/output) integer ip[natoms]
      , float   *const fp=nullptr // quantities (input)        float   fp[natoms or less]
      , double  *const *const dpp=nullptr // quantities (input/output) double* dpp[natoms]
    ) {
        auto const use = control::get("use.live.atom", 1.);
        if (0 == use) {
         // warn("single_atom::atom_update deactivated", 0); 
            return 0;
        } else {
            auto const stat = single_atom::atom_update(what, natoms, dp, ip, fp, dpp);
            if (stat) warn("single_atom::atom_update(%s, natoms=%d, ...) returned status= %i", what, natoms, int(stat));
            return stat;
        }
    } // live_atom_update

    inline double rho_Thomas_Fermi(double const Ekin, double & derivative) {
        if (Ekin <= 0) { derivative = 0; return 0; }
        double constexpr by_pi2 = 1./pow2(constants::pi);
        derivative = by_pi2 * std::sqrt(2*Ekin);
        double const value = derivative * Ekin * (2./3.);
        return value;
    } // rho_Thomas_Fermi


    status_t new_density_Thomas_Fermi(
          double rho_new[] // resulting density
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential
        , size_t const nall // number of all grid points
        , MPI_Comm const comm // grid communicator
        , double const n_electrons // required total number of electrons 
        , double const dV // grid volume element
        , int const echo // log level
        , int const maxiter=99 // maximum number of iterations
        , float const threshold=1e-8 // convergence criterion
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
    void block_average(real_t v4[4*4*4], real_t const v8[8*8*8], simple_stats::Stats<double> *s=nullptr) {
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


    double integrate_r2grid(double const r2coeff[], int32_t const nr2=4096, float const ar2=16.f) {
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
        for (size_t i = 0; i < all; ++i) {
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
        , float const r_cut=16.f // truncation radius of atom-cube interactions
    ) {
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

            // this should be the only loop over all atoms in the entire code, keep it free of I/O and only minimum number of operations
{ SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances of all atoms", echo);
            for (int32_t ia{0}; ia < n_all_atoms; ++ia) {
                auto const *const pos = xyzZ_all[ia];
                for (int ip{0}; ip < n_periodic_images; ++ip) {
                    auto const *const img = periodic_images[ip];
                    auto const dist2 = pow2(pos[0] + img[0] - cow[0])
                                    + pow2(pos[1] + img[1] - cow[1])
                                    + pow2(pos[2] + img[2] - cow[2]);
                    if (dist2 < r2search) {
                        double const img_pos[] = {pos[0] + img[0], pos[1] + img[1], pos[2] + img[2]};
                        double const Z = pos[3]; // pos points into xyzZ array
                        atom_images.push_back(atom_image_t(img_pos, Z, ia, periodic_shift[ip]));
                    } // inside search radius
                } // ip
            } // ia
} // timer
            if (echo > 2) std::printf("# rank#%i finds %lld atom images inside a %g %s search radius\n",
                                              me, atom_images.size(), r_search*Ang, _Ang);
        } // scope


        // now refine the search checking the proximity
        auto const r_close = r_cut + r_circum, r2close = pow2(r_close);

     // view2D<float> atom_block_distance2(atom_images.size(), n_blocks, 9e9);
        std::vector<uint32_t> n_close_blocks(atom_images.size(), 0),
                              n_close_atoms(n_blocks, 0);
        size_t distances_close{0}; // how many atom images are close to the domain

{ SimpleTimer timer(strip_path(__FILE__), __LINE__, "computing distances", echo);
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

        return atom_images;
    } // get_neigborhood


    std::vector<uint32_t> find_unique_atoms(std::vector<atom_image_t> & atom_images, int const echo=0) {
        // count how many periodic images each unique atom has
        std::map<int32_t,uint32_t> npi_map;
        for (auto const & ai : atom_images) {
            auto const global_atom_id = ai.atom_id_;
            ++npi_map[global_atom_id];
        } // ai
        // the size of this map determines how many different global_ids are involved
        auto const n_unique = npi_map.size();
        if (echo > 0) std::printf("# found %ld different atoms in %ld atom images\n", n_unique, atom_images.size());

        std::vector<uint32_t> global_atom_ids(n_unique); // prepare function result
        std::map<int32_t,uint32_t> id_map; // this map translates from global ids to 
        { // scope:
            simple_stats::Stats<float> s(0);
            uint32_t i_unique{0};
            for (auto const & ia : npi_map) {
                auto const global_atom_id = ia.first;
                id_map[global_atom_id] = i_unique;
                assert(global_atom_id >= 0);
                global_atom_ids[i_unique] = global_atom_id;
                ++i_unique;

                s.add(ia.second); // statistics about the number of periodic images
            } // ia
            assert(n_unique == i_unique);
            auto const s_interval = s.interval();
            if (echo > 3) std::printf("# distribution of periodic images is %s\n", s_interval.c_str());
        } // scope

        // translate the members atom_id_ in atom_images from global ids into unique ids
        for (auto ai : atom_images) {
            auto const global_atom_id = ai.atom_id_;
            auto const i_unique = id_map[global_atom_id];
            ai.atom_id_ = i_unique;
        } // ai

        return global_atom_ids;
    } // find_unique_atoms








    template<typename real_t>
    double add_r2grid_to_block(
          real_t values[8*8*8]
        , double const block_coords[3]
        , double const h[3] // Cartesian grid spacings
        , double const atom_center[3]
        , double const r2coeff[] // radial function on an r^2-grid
        , double const factor=1
        , int const echo=0
        , float const r_cut=-1 // -1: auto
        , float const ar2=16.f // r^2-grid parameter
        , int32_t const nr2=4096 // r^2-grid size
    ) {
        assert(ar2 > 0);
        assert(ar2 == 16.f);
        double const r_max = std::sqrt((nr2 - 1.)/ar2); // largest radius of the r^2-grid
        double const rcut = (-1 == r_cut) ? r_max : std::min(double(r_cut), r_max);
        double const r2cut = pow2(rcut);
        assert(ar2*r2cut < nr2);
        double added_charge{0};
        size_t grid_points_inside{0};
        if (1) { // ToDo: check if any of the radii can be inside
            for (int iz{0}; iz < 8; ++iz) {
            for (int iy{0}; iy < 8; ++iy) {
            for (int ix{0}; ix < 8; ++ix) {
                auto const r2 = pow2(block_coords[0] + h[0]*(ix + .5) - atom_center[0])
                              + pow2(block_coords[1] + h[1]*(iy + .5) - atom_center[1])
                              + pow2(block_coords[2] + h[2]*(iz + .5) - atom_center[2]);
                if (r2 < r2cut) {
                    int const ir2 = int(ar2*r2);

                    if (ir2 + 1 < nr2) {
                        double const w8 = ar2*r2 - ir2; // linear interpolation weight
                        auto const value_to_add = r2coeff[ir2] * (1 - w8) + r2coeff[ir2 + 1]*w8;
                        int const izyx = (iz*8 + iy)*8 + ix;
                        values[izyx] += factor*value_to_add;
                        added_charge += factor*value_to_add;
                    }
                    ++grid_points_inside;
                } // inside
            }}} // ix iy iz
        } // 1
        if (grid_points_inside && echo > 11) std::printf("# %lld grid points inside %g %s for grid block at [%g %g %g] Bohr\n",
            grid_points_inside, rcut*Ang, _Ang, block_coords[0], block_coords[1], block_coords[2]);
        return added_charge;
    } // add_r2grid_to_block



    template <int PROJECT0_OR_ADD1, int n8=8>
    size_t sho_project0_or_add1(
          double coeff[] // [nSHO(numax)] input if adding, result if projecting, coefficients are zyx-ordered
        , double values[] // grid block [n8*n8*n8], result if adding
        , double const block_coords[3] // block origin
        , double const h[3] // grid spacing
        , double const atom_center[3] // where is the atom image
        , int const numax // how many
        , double const sigma // Gaussian spread
        , float const r_cut
        , int const echo=0 // log-level
    ) {
        // auto const rcut = sho_projection::truncation_radius(sigma, numax);
        float const r2cut = pow2(r_cut);
        assert(sigma > 0);
        auto const sigma_inv = 1./sigma;
        int const nSHO = sho_tools::nSHO(numax);

        // ToDo: analyze if the grid spacing is small enough for this \sigma

        int const n1HO = sho_tools::n1HO(numax);
        view3D<double> H1d(3,n8,n1HO);
        float radius2[3][n8];
        for (int dir = 0; dir < 3; ++dir) {
            if (echo > 55) std::printf("\n# Hermite polynomials for %c-direction:\n", 'x' + dir);
            for (int i8 = 0; i8 < n8; ++i8) {
                auto const x = block_coords[dir] + (i8 + .5)*h[dir] - atom_center[dir];
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

#ifdef    DEVEL
        if (1 == PROJECT0_OR_ADD1) {
            if (echo > 16) {
                std::printf("# addition coefficients ");
                for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                    std::printf(" %g", std::real(coeff[iSHO]));
                } // iSHO
                std::printf("\n");
            } // echo
        } // ADD
#endif // DEVEL


        size_t nhits{0};
        for (int iz = 0; iz < n8; ++iz) {
        for (int iy = 0; iy < n8; ++iy) {
        for (int ix = 0; ix < n8; ++ix) {
            int const ixyz = (iz*n8 + iy)*n8 + ix;

            double val{0};
            if (radius2[0][ix] + radius2[1][iy] + radius2[2][iz] < r2cut) {
//                    if (echo > 6) std::printf("%g %g\n", std::sqrt(vz*vz + vy*vy + vx*vx), val); // plot function value vs r
                if (0 == PROJECT0_OR_ADD1) {
                    val = values[ixyz]; // load
                } // project
                int iSHO{0};
                for (int nz = 0; nz <= numax; ++nz) {                    auto const H1d_z = H1d(2,iz,nz);
                    for (int ny = 0; ny <= numax - nz; ++ny) {           auto const H1d_y = H1d(1,iy,ny);
                        for (int nx = 0; nx <= numax - nz - ny; ++nx) {  auto const H1d_x = H1d(0,ix,nx);
                            auto const H3d = H1d_z * H1d_y * H1d_x;
                            if (1 == PROJECT0_OR_ADD1) {
                                val += coeff[iSHO] * H3d; // here, the addition happens
                            } else {
                                coeff[iSHO] += val * H3d; // here, the projection happens
                            }
                            ++iSHO; // in sho_tools::zyx_order
                        } // nx
                    } // ny
                } // nz
                assert( nSHO == iSHO );
                ++nhits;
            } // inside radius

            if (1 == PROJECT0_OR_ADD1) {
                values[ixyz] += val; // load-modify-store, must be atomic if threads are involved
            } // write back (add)

        }}} // ix iy iz

#ifdef    DEVEL
        if (0 == PROJECT0_OR_ADD1) {
            if (echo > 16) {
                std::printf("# projection coefficients ");
                for (int iSHO = 0; iSHO < nSHO; ++iSHO) {
                    std::printf(" %g", std::real(coeff[iSHO]));
                } // iSHO
                std::printf("\n");
            } // echo
        } // PROJECT
#endif // DEVEL

        return nhits;
    } // sho_project0_or_add1


    template <int n8=8>
    size_t sho_project(
          double coeff[] // result [nSHO(numax)], coefficients are zyx-ordered
        , double const values[] // grid block [n8*n8*n8], result if adding
        , double const block_coords[3] // block origin
        , double const h[3] // grid spacing
        , double const center[3] // where is the atom image
        , int const numax // how many
        , double const sigma // Gaussian spread
        , float const r_cut
        , int const echo=0 // log-level
    ) {
        auto constexpr PROJECT0_OR_ADD1 = 0;
        return sho_project0_or_add1<PROJECT0_OR_ADD1,n8>(coeff, (double*)values,
                          block_coords, h, center, numax, sigma, r_cut, echo);
    } // sho_project

    template <int n8=8>
    size_t sho_add(
          double values[] // result grid block [n8*n8*n8] to add to
        , double const coeff[] // input [nSHO(numax)], coefficients are zyx-ordered
        , double const block_coords[3] // block origin
        , double const h[3] // grid spacing
        , double const center[3] // where is the atom image
        , int const numax // how many
        , double const sigma // Gaussian spread
        , float const r_cut
        , int const echo=0 // log-level
    ) {
        auto constexpr PROJECT0_OR_ADD1 = 1;
        return sho_project0_or_add1<PROJECT0_OR_ADD1,n8>((double*)coeff, values,
                          block_coords, h, center, numax, sigma, r_cut, echo);
    } // sho_add

    status_t add_to_grid(
          view2D<double> & grid_values // result: add to these grid array values, usually augmented_density
        , view2D<double> const & block_coords
        , uint32_t const n_blocks
        , data_list<double> const & atom_coeff // [natoms][nSHO]
        , std::vector<int32_t> const & lmax // lmax_qlm[natoms]
        , std::vector<double> const & sigma // sigma_cmp[natoms]
        , std::vector<atom_image_t> const & atom_images
        , double const grid_spacings[3] // h
        , int const echo=0
    ) {
        std::vector<size_t> hits_per_atom(sigma.size(), 0);
        for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // OMP PARALLEL
            size_t hits_per_block{0};
            for (auto const & ai : atom_images) {
                auto const iatom = ai.atom_id_;
                float const r_cut = 9*sigma[iatom];
                // auto const hits = sho_add(grid_values[ilb], atom_coeff[iatom], block_coords[ilb],
                //               grid_spacings, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                auto const hits = sho_project0_or_add1<1,8>((double*)atom_coeff[iatom], grid_values[ilb], block_coords[ilb],
                              grid_spacings, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                hits_per_atom[iatom] += hits;
                hits_per_block       += hits;
            } // ai
            if (echo > 11 && hits_per_block > 0) std::printf("# %s: block #%i at %g %g %g Bohr has %.3f k hits\n", __func__,
                ilb, block_coords[ilb][0], block_coords[ilb][1], block_coords[ilb][2], hits_per_block*1e-3);
        } // ilb
        for (size_t iatom{0}; iatom < sigma.size()*(echo > 0); ++iatom) {
            if (hits_per_atom[iatom] > 0) std::printf("# %s: atom #%i has %.3f k hits\n", __func__, iatom, hits_per_atom[iatom]*1e-3);
        } // iatom
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
        , double const grid_spacings[3] // h
        , int const echo=0
    ) {
        auto const natoms = sigma.size();
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
                // auto const hits = sho_project(atom_coeff[iatom], grid_values[ilb], block_coords[ilb],
                //               grid_spacings, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                auto const hits = sho_project0_or_add1<0,8>(atom_coeff[iatom], (double*)grid_values[ilb], block_coords[ilb],
                              grid_spacings, ai.pos_, lmax[iatom], sigma[iatom], r_cut, echo);
                hits_per_block[ilb] += hits;
                hits_per_atom       += hits;
            } // ilb
            if (echo > 0 && hits_per_atom > 0) std::printf("# %s: image of atom #%i at %g %g %g Bohr has %.3f k hits\n", __func__, 
                iatom, ai.pos_[0], ai.pos_[1], ai.pos_[2], hits_per_atom*1e-3);
        } // ai
        for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) {
            if (echo > 11 && hits_per_block[ilb] > 0) std::printf("# %s: block #%i at %g %g %g Bohr has %.3f k hits\n", __func__,
                ilb, block_coords[ilb][0], block_coords[ilb][1], block_coords[ilb][2], hits_per_block[ilb]*1e-3);
        } // ilb
        return 0;
    } // project_grid



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
    ) {
        auto const me = mpi_parallel::rank(comm);
        // ToDo: get atom_r2coeff from atom owners

        double rho_total{0}; // prepare result
        std::vector<double> added_charge_per_atom(natoms, 0.0); // DEBUG
        for (auto const & ai : atom_images) {
            auto const iatom = ai.atom_id_;
            assert(0 <= iatom); assert(iatom < natoms);
            if (echo > 12) std::printf("# rank#%i adds to %s for atom#%i at position %g %g %g %s\n", 
                    me, quantity_name, iatom, ai.pos_[0]*Ang, ai.pos_[1]*Ang, ai.pos_[2]*Ang, _Ang);
            double added_charge{0};
            for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // local blocks
                added_charge += add_r2grid_to_block(grid_quantity[ilb], block_coords[ilb], g.grid_spacings(), ai.pos_, atom_r2coeff[iatom], factor, echo);
            } // ilb
            if (echo > 13) std::printf("# rank#%i atom#%i %g electrons added for periodic image %i %i %i\n",
                                me, iatom, added_charge*g.dV(), ai.shifts_[0], ai.shifts_[1], ai.shifts_[2]);
            rho_total += added_charge;
            added_charge_per_atom[iatom] += added_charge; // DEBUG
        } // ia
        rho_total = mpi_parallel::sum(rho_total, comm)*g.dV();
        if (echo > 9) { // DEBUG
            double rho_total_checksum{0};
            for (int iatom{0}; iatom < natoms; ++iatom) {
                std::printf("# rank#%i atom#%i  %g electrons added\n",
                      me, iatom, added_charge_per_atom[iatom]*g.dV());
                rho_total_checksum += added_charge_per_atom[iatom];
            } // iatom
            std::printf("# rank#%i added %g electrons to %s\n", me, rho_total_checksum*g.dV(), quantity_name);
        } // echo // DEBUG
        if (echo > 3) std::printf("# added %g electrons to %s\n", rho_total, quantity_name);
        return rho_total;
    } // add_r2grid_quantity




    status_t SCF(int const echo=0) {
        status_t stat(0);

        auto const comm = MPI_COMM_WORLD;
        auto const me = mpi_parallel::rank(comm);
        auto const nprocs = mpi_parallel::size(comm);

        // load geometry from +geometry.file=atoms.xyz (a good candidate is +geometry.file=test/geo/Cu40Zr22.xyz)
        real_space::grid_t g; // entire grid descriptor
        view2D<double> xyzZ_all;  // coordinates for all atoms
        int32_t n_all_atoms;  // number of all atoms

        if (0 == me) { // MPI master task
            auto const stat_init = self_consistency::init_geometry_and_grid(g, xyzZ_all, n_all_atoms, 8, echo);
            if (stat_init) warn("init_geometry_and_grid returned status= %i", int(stat_init));
            stat += stat_init;
        } // master
        mpi_parallel::broadcast(&n_all_atoms, comm);
        if (0 != me) { xyzZ_all = view2D<double>(n_all_atoms, 4); }
        mpi_parallel::broadcast(xyzZ_all.data(), comm, n_all_atoms*4);
        mpi_parallel::broadcast((char*)&g, comm, sizeof(g));

        // create a coarse grid descriptor
        real_space::grid_t gc(g[0]/2, g[1]/2, g[2]/2);
        {
            gc.set_boundary_conditions(g.boundary_conditions());
            gc.set_cell_shape(g.cell, echo*0); // muted
            gc.set_grid_spacing(g.h[0]*2, g.h[1]*2, g.h[2]*2);
            for (int d{0}; d < 3; ++d) { assert(0 == (gc[d] & 0x3)); } // all grid numbers must be a multiple of 4
            auto const max_grid_spacing = std::max(std::max(std::max(1e-9, gc.h[0]), gc.h[1]), gc.h[2]);
            if (echo > 1) std::printf("# use  %g %g %g  %s coarse grid spacing, corresponds to %.1f Ry\n",
                      gc.h[0]*Ang, gc.h[1]*Ang, gc.h[2]*Ang, _Ang, pow2(constants::pi/max_grid_spacing));
        }
        double const grid_center[] = {g[0]*g.h[0]*.5, g[1]*g.h[1]*.5, g[2]*g.h[2]*.5}; // reference point for atomic positions

        // distribute the dense grid in 8x8x8 grid blocks to parallel owners
        parallel_poisson::parallel_grid_t pg(g, comm, 8, echo);
        if (echo > 1) { auto const nb = pg.grid_blocks(); std::printf("# use  %d %d %d  grid blocks\n", nb[0], nb[1], nb[2]); }

        auto const n_blocks = pg.n_local();
        view2D<double> block_coords(n_blocks, 4, 0.0);
        auto atom_images = get_neighborhood(block_coords, n_all_atoms, xyzZ_all, pg, g, grid_center, echo);
        // for now atom_images.atom_id_ are in [0, n_all_atoms)

        auto const global_atom_ids = find_unique_atoms(atom_images, echo);
        // from here atom_images.atom_id_ are in [0, global_atom_ids.size())



        // distribute the atom ownership TODO
        // simple distribution model: owner_rank == atom_id % nprocs, advantage: every rank can compute the owner
        // alternativ: make a weight grid pg.nb^3, init 0, add weight where atoms are, call load_balancer



        // int32_t const na = n_all_atoms; // number of owned atoms          TODO atom ownership not parallelized TODO
        int32_t const na = global_atom_ids.size(); assert(global_atom_ids.size() == na);
        uint32_t const natoms = global_atom_ids.size(); // number of locally relevant atoms

        float take_atomic_valence_densities{1}; // 100% of the smooth spherical atomic valence densities is included in the smooth core densities

        std::vector<double> sigma_cmp(na, 1.); // spread of the Gaussian used in the compensation charges
        char const pawdata_from = (*control::get("pawdata.from", "auto")) | 32; // 'a': auto generate, 'f': pawxml_import
        std::vector<int32_t> numax(na, ('f' == pawdata_from)?-9:-1); // -1: LivePAW, -9: load from pawxml files
        std::vector<int32_t> lmax_qlm(na, -1);
        std::vector<int32_t> lmax_vlm(na, -1);

        // initialize and get sigma, lmax for each atom
        std::vector<int32_t> n_atom_rho(na, 0);
        data_list<double> atom_qlm, atom_vlm, atom_rho, atom_mat;
        data_list<double> atoms_qlm, atoms_vlm;
        {
            std::vector<double> Za(na, 0.);
            for (int32_t ia{0}; ia < na; ++ia) { Za[ia] = xyzZ_all(global_atom_ids[ia],3); }
            std::vector<float> ionization(na, 0.f);

            stat += live_atom_update("initialize", na, Za.data(), numax.data(), ionization.data(), (double**)1);
            stat += live_atom_update("lmax qlm",   na,    nullptr, lmax_qlm.data(), &take_atomic_valence_densities);
            stat += live_atom_update("lmax vlm",   na, (double*)1, lmax_vlm.data());
            stat += live_atom_update("sigma cmp",  na, sigma_cmp.data());

            view2D<uint32_t> num(4, na, 0); // how many
            for (int32_t ia{0}; ia < na; ++ia) {
                num(0,ia) = pow2(1 + lmax_qlm[ia]); // qlm
                num(1,ia) = pow2(1 + lmax_vlm[ia]); // vlm
                int const ncoeff = sho_tools::nSHO(numax[ia]);
                num(2,ia) =   pow2(ncoeff);         // aDm
                num(3,ia) = 2*pow2(ncoeff);         // aHm+aSm
                n_atom_rho[ia] = num(2,ia); // number of coefficients for the atomic density matrices
            } // ia
            // get memory in the form of data_list containers
            atom_qlm = data_list<double>(na, num[0], 0.0);
            atom_vlm = data_list<double>(na, num[1], 0.0);
            atom_rho = data_list<double>(na, num[2], 0.0);
            atom_mat = data_list<double>(na, num[3], 0.0);

            // for addition and projection
            view2D<uint32_t> numb(2, natoms, 0); // how many
            for (int32_t iatom{0}; iatom < natoms; ++iatom) {
                numb(0,iatom) = sho_tools::nSHO(lmax_qlm[iatom]);
                numb(1,iatom) = sho_tools::nSHO(lmax_vlm[iatom]);
            } // iatom
            atoms_qlm = data_list<double>(natoms, numb[0], 0.0);
            atoms_vlm = data_list<double>(natoms, numb[1], 0.0);

        } // scope
        xyzZ_all = view2D<double>(0, 0, 0.0); // clear, xyzZ_all should not be used after this

        // the r^2-grids are used to bring radial quantities to a Cartesian grid
        std::vector<int32_t> nr2(na, 1 << 12); // 4096
        data_list<double> atom_vbar(nr2, 0.0); // zero potentials
        data_list<double> atom_rhoc(nr2, 0.0); // core densities



        // allocate CPU memory for grid arrays
        view2D<double> augmented_density(n_blocks, 8*8*8, 0.0);
        view2D<double> V_electrostatic(  n_blocks, 8*8*8, 0.0);
        view2D<double> core_density(     n_blocks, 8*8*8, 0.0);
        view3D<double> valence_density(2,n_blocks, 8*8*8, 0.0); // 0: density, 1: energy derivative w.r.t. the Fermi level


        double constexpr Y00 = .28209479177387817; // == 1/sqrt(4*pi)
        double constexpr Y00sq = pow2(Y00);        // == 1/(4*pi)

        { // scope: collect start densities from atomic densities

            // get smooth core densities on r^2-grids
            stat += live_atom_update("core densities", na, 0, nr2.data(), 0, atom_rhoc.data());

            if (take_atomic_valence_densities > 0) {
                auto & atom_rhov = atom_vbar; // use memory of vbar for the moment
                stat += live_atom_update("valence densities", na, 0, nr2.data(), 0, atom_rhov.data());
                for (int32_t ia{0}; ia < na; ++ia) {
                    // add valence density to r^2-gridded core density
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core electrons\n",         me, ia, integrate_r2grid(atom_rhoc[ia], nr2[ia]));
                    add_product(atom_rhoc[ia], nr2[ia], atom_rhov[ia], take_atomic_valence_densities*1.);
                    if (echo > 15) std::printf("# rank#%i atom#%i wants to add %g core+valence electrons\n", me, ia, integrate_r2grid(atom_rhoc[ia], nr2[ia]));
                } // ia
            } // take_atomic_valence_densities > 0

        } // scope

        // These data items are needed for the Kohn-Sham Hamiltonian
        // stat += live_atom_update("projectors", na, sigma_prj.data(), numax.data());

        // configure the Poisson solver
        auto  const es_method    = control::get("poisson.method", "cg"); // {cg, sd}
        int   const es_echo      = control::get("poisson.echo", echo*1.);
        float const es_threshold = control::get("poisson.threshold", 3e-8);
        int   const es_maxiter   = control::get("poisson.maxiter", 200.);
        int   const es_miniter   = control::get("poisson.miniter", 3.);
        int   const es_restart   = control::get("poisson.restart", 4096.);

        // configure the self-consistency loop
        int   const scf_maxiter   = control::get("scf.maxiter", 1.);

        double E_Fermi{0}; // Fermi level

        double nve{0}; // non-const total number of valence electrons
        { // scope: determine the number of valence electrons
            auto const keyword_valence_electrons = "valence.electrons";
            if ('a' == (*control::get(keyword_valence_electrons, "auto") | 32)) {
                std::vector<double> n_electrons_a(na, 0.); // number of valence electrons added by each atom
                stat += live_atom_update("#valence electrons", na, n_electrons_a.data());
                nve = std::accumulate(n_electrons_a.begin(), n_electrons_a.end(), 0.0);
                if (echo > 0) std::printf("\n# %s=auto --> %g valence electrons\n\n", keyword_valence_electrons, nve);
                control::set(keyword_valence_electrons, nve, control::echo_set_without_warning);
            } else {
                nve = control::get(keyword_valence_electrons, 0.0);
                if (echo > 0) std::printf("\n# %s=%g\n\n", keyword_valence_electrons, nve);
            } // auto
        } // scope
        double const n_valence_electrons = nve; // do not use nve beyond this point

        sho_unitary::Unitary_SHO_Transform const unitary(9);

        int scf_iteration{0};
        bool scf_run{true};
        do { // self-consistency loop

            mpi_parallel::barrier(comm);
            ++scf_iteration;
            if (echo > 3) std::printf("#\n# Start SCF iteration #%i\n#\n", scf_iteration);


            auto const total_charge_added = add_r2grid_quantity(core_density, "smooth core density", atom_rhoc,
                                            atom_images, natoms, block_coords, n_blocks, g, comm, echo, Y00sq);
            if (echo > 0) std::printf("# %g electrons added as smooth core density\n", total_charge_added);

            // compose density
            set(augmented_density[0], n_blocks*size_t(8*8*8), core_density[0]);
            if (take_atomic_valence_densities < 1) { // less than 100% atomic valence densities
                add_product(augmented_density[0], n_blocks*size_t(8*8*8), valence_density(0,0), 1.);
            } // take_atomic_valence_densities < 1

            print_stats(augmented_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# smooth density");

            view2D<double> V_xc(n_blocks, 8*8*8, 0.0);
            { // scope: eval the XC potential and energy
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
                if (echo > 2) std::printf("# exchange-correlation energy on grid %.9f %s, double counting %.9f %s\n", E_xc*eV,_eV, E_dc*eV,_eV);
                // grid_xc_energy = E_xc;
            } // scope
            print_stats(V_xc[0], n_blocks*size_t(8*8*8), comm, echo > 1, 0, "# smooth exchange-correlation potential", eV, _eV);

            stat += live_atom_update("qlm charges", na, 0, 0, 0, atom_qlm.data());

            for (int iatom{0}; iatom < natoms; ++iatom) { // ToDo: atom-owner broadcasts denormalized qlms
                auto const ia = iatom; assert(na == natoms);
                auto const stat_den = sho_projection::denormalize_electrostatics(atoms_qlm[iatom], atom_qlm[ia], lmax_qlm[ia], sigma_cmp[ia], unitary, echo);
                if (stat_den) warn("denormalize_electrostatics failed with status= %i", int(stat_den));
                stat += stat_den;
            } // iatom
            add_to_grid(augmented_density, block_coords, n_blocks, atoms_qlm, lmax_qlm, sigma_cmp, atom_images, g.grid_spacings(), echo);

            print_stats(augmented_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# smooth augmented_density");


            { // scope: Poisson equation
                if (echo > 3) std::printf("#\n# Solve the Poisson equation iteratively with %d ranks in SCF iteration #%i\n#\n", nprocs, scf_iteration);
                float es_residual{0};
                auto const es_stat = parallel_poisson::solve(V_electrostatic[0], augmented_density[0],
                                                pg, *es_method, es_echo, es_threshold,
                                                &es_residual, es_maxiter, es_miniter, es_restart);
                if (echo > 2) std::printf("# Poisson equation %s, residual= %.2e a.u.\n#\n", es_stat?"failed":"converged", es_residual);
                stat += es_stat;
                if (es_stat) warn("parallel_poisson::solve returned status= %i", int(es_stat));
            } // scope

            print_stats(V_electrostatic[0], n_blocks*size_t(8*8*8), comm, echo > 0, 0, "# smooth electrostatic potential", eV, _eV);

            // project the electrostatic grid onto the localized compensation charges, TODO: multiply with g.dV()
            project_grid(atoms_vlm, V_electrostatic, block_coords, n_blocks, lmax_vlm, sigma_cmp, atom_images, g.grid_spacings(), echo);
            // ToDo: allreduce atoms_vlm and normalize to lm-representation at atom-owner

            for (int ia{0}; ia < na; ++ia) { // ToDo: addreduce vlm in atom-owner process
                auto const iatom = ia; assert(na == natoms);
                scale(atoms_vlm[iatom], sho_tools::nSHO(lmax_vlm[ia]), g.dV()); // scaling with g.dV() is NOT performed inside project_grid()
                auto const stat_ren = sho_projection::renormalize_electrostatics(atom_vlm[ia], atoms_vlm[iatom], lmax_vlm[ia], sigma_cmp[ia], unitary, echo);
                if (stat_ren) warn("renormalize_electrostatics failed with status= %i", int(stat_ren));
                stat += stat_ren;
            } // iatom

#ifdef    DEVEL
            if (echo > 3) {
                for (int ia{0}; ia < na; ++ia) {
                    std::printf("# potential projection for atom #%d v_00 = %.9f %s\n", ia, atom_vlm[ia][00]*Y00*eV, _eV);
                    int const ellmax_show = std::min(lmax_vlm[ia], 2);
                    for (int ell = 1; ell <= ellmax_show; ++ell) {
                        double const unitfactor = Y00 * eV * std::pow(Ang, -ell);
                        int const ilm0 = sho_tools::lm_index(ell, -ell);
                        std::printf("# potential projection for atom #%d v_%im =", ia, ell);
                        printf_vector(" %.6f", &atom_vlm[ia][ilm0], 2*ell + 1, nullptr, unitfactor);
                        std::printf(" %s %s^%i\n", _eV, _Ang, -ell);
                    } // ell
                } // ia
            } // echo
#endif // DEVEL

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

            add_r2grid_quantity(V_effective, "smooth effective potential", atom_vbar,
                                atom_images, natoms, block_coords, n_blocks, g, comm, echo*0, Y00);

            print_stats(V_effective[0], n_blocks*size_t(8*8*8), comm, echo > 0, 0, "# smooth effective potential", eV, _eV);

            view2D<double> V_coarse(n_blocks, 4*4*4, 0.0);
            { // scope: reduce them to 4x4x4 grid points per block
                for (uint32_t ilb{0}; ilb < n_blocks; ++ilb) { // parallel loop over local blocks
                    block_average(V_coarse[ilb], V_effective[ilb]);
                } // ilb
            } // scope
            print_stats(V_coarse[0], n_blocks*size_t(4*4*4), comm, echo > 0, 0, "# coarse effective potential", eV, _eV);



            // ToDo: call energy-contour integration to find a new density
            view2D<double> new_valence_density(n_blocks, 8*8*8, 0.0);
            { // scope: apply Thomas-Fermi approximation
                auto const stat_TF = new_density_Thomas_Fermi(new_valence_density[0], E_Fermi, V_effective[0],
                                n_blocks*size_t(512), comm, n_valence_electrons, g.dV(), echo);
                stat += stat_TF;
                if (stat_TF) warn("# new_density_Thomas_Fermi returned status= %i", int(stat_TF));
                print_stats(new_valence_density[0], n_blocks*size_t(8*8*8), comm, echo > 0, g.dV(), "# new Thomas-Fermi density");
            } // scope



            float rho_mixing_ratios[] = {.5, .5, .5}; // for spherical {core, semicore, valence} density
            stat += live_atom_update("atomic density matrices", na, 0, 0, rho_mixing_ratios, atom_rho.data());

            scf_run = (scf_iteration < scf_maxiter);

        } while (scf_run); // self-consistency loop 

        stat += live_atom_update("memory cleanup", na);
        return stat;
    } // SCF

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
        stat += test_scf(echo);
        if (!already_initialized) mpi_parallel::finalize();
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace parallel_potential
