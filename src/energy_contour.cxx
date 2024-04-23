// This file is part of AngstromCube under MIT License

#include <cstdint> // int8_t
#include <cassert> // assert
#include <cstdio> // std::printf, ::snprintf
#include <vector> // std::vector<T>
#include <complex> // std::complex<real_t>
#include <algorithm> // std::min, ::max

#include "energy_contour.hxx"

#include "control.hxx" // ::get
#include "display_units.h" // eV, _eV, Ang, _Ang, Kelvin, _Kelvin
#include "mpi_parallel.hxx" // ::rank, ::comm, MPI_COMM_WORLD
#include "data_view.hxx" // view2D<T>
#include "parallel_poisson.hxx" // ::parallel_grid_t
#include "action_plan.hxx" // action_plan_t
#include "green_function.hxx" // ::construct_Green_function, ::update_atom_matrices, ::update_phases, ::update_energy_parameter, ::update_potential
#include "real_space.hxx" // ::grid_t
#include "data_list.hxx" // data_list<T>
#include "recorded_warnings.hxx" // warn
#include "inline_math.hxx" // set, add_product
#include "green_solver.hxx" // green_solver_t
#include "recorded_warnings.hxx" // error, warn
#include "sho_tools.hxx" // ::nSHO
#include "sho_projection.hxx" // ::get_sho_prefactors
#include "brillouin_zone.hxx" // ::get_kpoint_mesh, ::WEIGHT

namespace energy_contour {

    Integrator::Integrator( // implementation of constructor
          real_space::grid_t const & gc // coarse grid descriptor
        , std::vector<double> const & xyzZinso // all atoms
        , int const echo // verbosity
    ) {
        if (echo > 0) std::printf("# construct %s with grid=[%d %d %d]\n", __func__, gc[0], gc[1], gc[2]);
        plan_ = new action_plan_t(); // CPU memory for the plan
        auto const stat = green_function::construct_Green_function(*plan_,
                            gc.grid_points(), gc.boundary_conditions(), gc.grid_spacings(),
                            xyzZinso, echo);
        if (stat) warn("construct_Green_function returned status= %i", int(stat));

        if (echo > 0) std::printf("# move green_solver_t\n");
        solver_ = new green_solver_t(plan_, echo);
        if (echo > 0) std::printf("# constructed %s\n", __func__);
    } // constructor

    typedef std::complex<double> Complex;

    status_t show_contour(std::vector<Complex> const & E_points, int const echo=0) {
        // an ASCII plot where the energy points are located
        if (echo < 1) return 0;
        auto const nE = E_points.size();
        if (nE < 1) return 0;
        float ex[2][2] = {{9e9f, -9e9f}, {9e9f, -9e9f}};
        for (auto const ep : E_points) {
            ex[0][0] = std::min(ex[0][0], float(ep.real()));
            ex[0][1] = std::max(ex[0][1], float(ep.real()));
            ex[1][0] = std::min(ex[1][0], float(ep.imag()));
            ex[1][1] = std::max(ex[1][1], float(ep.imag()));
        } // ep
        std::printf("#\n# energy contour %ld points within [%g, %g %s] and [%g, %g %s]\n",
            nE, ex[0][0]*eV, ex[0][1]*eV, _eV, ex[1][0]*Kelvin, ex[1][1]*Kelvin, _Kelvin);
        return 0;
    } // show_contour

    std::vector<Complex> get_energy_mesh(
          std::vector<Complex> & w8
        , double   const kBT
        , double   const eBot
        , unsigned const nBot
        , unsigned const nPar
        , unsigned const nFer
        , int      const nPol
        , int const echo=0
    ) {
        auto const nE = nBot + nPar + nFer + nPol;
        if (echo > 3) std::printf("# energy_contour with %d + %d + %d + %d poles = %d points\n", nBot, nPar, nFer, nPol, nE); 
        // from juKKR: https://iffgit.fz-juelich.de/kkr/jukkr/-/blob/master/source/KKRhost/emesht.f90

        assert(nPol > 0 && "Only nPol > 0 implemented");
        assert(kBT > 0);

        std::vector<Complex> Ep(0);
        w8.resize(0);
        if (nE > 0) {
            Ep.reserve(nE);
            w8.reserve(nE);
        }

        double constexpr E_mu = 0; // reference Fermi level is at zero
        auto const pikBT = constants::pi * kBT;
        double xi[128], wi[128];

        int jE{0};
        { // scope: sample a path from the real-valued band bottom to a point perpendicular to the real axis
            auto const *const path = "bottom";
            int const mBot = Gauss_Legendre_quadrature(xi, wi, nBot, echo);
            double const dE = nPol*pikBT;
            assert(dE > 0);
            for (int iE{0}; iE < mBot; ++iE) {
                Ep.push_back(Complex(eBot, (xi[iE] + 1.)*dE));
                w8.push_back(Complex(0.0, wi[iE]*dE)); // imaginary weights
                if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) %s, weight=(%g, %g)\n",
                    jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin, path, w8[jE].real(), w8[jE].imag());
                ++jE;
            } // iE
        } // scope

        { // scope: sample a path parallel to the real axis up to the Fermi level
            auto const *const path = "parallel";
            int const mPar = Gauss_Legendre_quadrature(xi, wi, nPar, echo);
            double const dE = (E_mu - 30*kBT - eBot)*0.5;
            for (int iE{0}; iE < mPar; ++iE) {
                Ep.push_back(Complex(eBot + (xi[iE] + 1.)*dE, nPol*2*pikBT));
                w8.push_back(Complex(wi[iE]*dE, 0.0)); // real-valued weights
                if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) %s, weight=(%g, %g)\n",
                    jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin, path, w8[jE].real(), w8[jE].imag());
                ++jE;
            } // iE
        } // scope

        { // scope: sample the Fermi-Dirac distribution from the Fermi level upwards in energy parallel to the real axis
            auto const *const path = "fermidirac";
            int const mFer = Gauss_Fermi_Dirac_quadrature(xi, wi, nFer, echo);
            double const dE = 30*kBT;
            for (int iE{0}; iE < mFer; ++iE) {
                Ep.push_back(Complex(E_mu + xi[iE]*dE, nPol*2*pikBT));
                w8.push_back(Complex(wi[iE]*dE, 0.0)); // real-valued weights
                if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) %s, weight=(%g, %g)\n",
                    jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin, path, w8[jE].real(), w8[jE].imag());
                ++jE;
            } // iE
        } // scope

        { // scope: add Matsubara poles
            auto const *const path = "matsubara";
            for (int iE{nPol}; iE > 0; --iE) {
                Ep.push_back(Complex(E_mu, (2*iE - 1)*pikBT));
                w8.push_back(Complex(0.0, -2*pikBT)); // imaginary weights
                if (echo > 8) std::printf("# energy_mesh[%2i]=(%11.6f %s, %6.1f %s) %s, weight=(%g, %g)\n",
                    jE, Ep[jE].real()*eV, _eV, Ep[jE].imag()*Kelvin, _Kelvin, path, w8[jE].real(), w8[jE].imag());
                ++jE;
            } // iE
        } // scope

        if (echo > 6) std::printf("# energy_mesh %2i points (expected %d points)\n", jE, nE);
        assert(Ep.size() == jE); assert(w8.size() == jE);

        show_contour(Ep, echo);
        return Ep;
    } // get_energy_mesh

    std::vector<Complex> get_energy_mesh(std::vector<Complex> & w8, int const echo=0) {
        double   const kBT  = control::get("energy_contour.temperature", 1e-2);
        double   const eBot = control::get("energy_contour.band.bottom", -1.); // w.r.t. the Fermi Level at 0
        unsigned const nBot = control::get("energy_contour.bottom",    5.);
        unsigned const nPar = control::get("energy_contour.parallel", 15.);
        unsigned const nFer = control::get("energy_contour.fermidirac", 3.);
        int      const nPol = control::get("energy_contour.matsubara", 3.);
        return get_energy_mesh(w8, kBT, eBot, nBot, nPar, nFer, nPol, echo);
    } // get_energy_mesh

    template <typename real_t>
    real_t sum(real_t const rho[], size_t const n) {
        real_t s{0};
        for (size_t i{0}; i < n; ++i) {
            s += rho[i];
        } // i
        return s;
    } // sum

    template <typename real_t>
    double maxval(real_t const rho[], size_t const n) {
        double s{0};
        for (size_t i{0}; i < n; ++i) {
            s = std::max(s, std::real(rho[i]));
        } // i
        return s;
    } // maxval

    status_t Integrator::integrate(
          double rho_888[] // resulting density in [nblocks][8*8*8] data layout
        , double & Fermi_level // Fermi level
        , double const Vtot[] // input potential in [nblocks][4*4*4]
        , data_list<double> const & atom_mat // atomic_Hamiltonian elements, only in atom owner ranks
        , std::vector<int32_t> const & numax_prj
        , std::vector<double> const & sigma_prj
        , parallel_poisson::parallel_grid_t const & pg // nblocks == pg.n_local()
        , double const n_electrons // =1 // required total number of electrons 
        , double const dV // =1 // grid volume element
        , int const echo // =0 // log level
    ) {
        status_t stat(0);
        size_t constexpr n4x4x4 = 4*4*4;
        auto const comm = mpi_parallel::comm(); // == MPI_COMM_WORLD
        auto const me = mpi_parallel::rank(comm);

        int const check = control::get("check", 0.);
        int const iterations = control::get("green_solver.iterations", 99.);
        if (echo > 0) std::printf("\n# energy_contour::integration(E_Fermi=%g %s, %g electrons, echo=%d) +check=%i\n", Fermi_level*eV, _eV, n_electrons, echo, check);

        auto const nblocks = pg.n_local();
        assert(nullptr != plan_);
        assert(nullptr != solver_);
        auto & plan = *plan_;
        plan.echo = echo >> 2; // lower verbosity

        if (plan.nCols != nblocks) warn("model assumes that each local block has one RHS, found n_local= %d and p.nRHS= %d", nblocks, plan.nCols);

        int constexpr Noco = 1;

        auto const nAtoms = atom_mat.nrows();
        std::vector<std::vector<double>> AtomMatrices(nAtoms);
        for (int iAtom{0}; iAtom < nAtoms; ++iAtom) {
            auto const nc2 = atom_mat.ncols(iAtom);
            auto const numax = numax_prj.at(iAtom);
            auto const nc = sho_tools::nSHO(numax);
            assert(nc2 == 2*nc*nc);
            auto const rescale = sho_projection::get_sho_prefactors(numax, sigma_prj.at(iAtom));
            AtomMatrices[iAtom] = std::vector<double>(nc2);
            for (int i{0}; i < nc; ++i) {
                auto const rescale_i = rescale.at(i);
                for (int j{0}; j < nc; ++j) {
                    // atom matrices need to be prepared for projection with unnormalized Gauss-Hermite functions
                    AtomMatrices[iAtom][(0*nc + i)*nc + j] = rescale_i * atom_mat[iAtom][(0*nc + i)*nc + j] * rescale[j]; // hmt
                    AtomMatrices[iAtom][(1*nc + i)*nc + j] = rescale_i * atom_mat[iAtom][(1*nc + i)*nc + j] * rescale[j]; // ovl
                } // j
            } // i
            if (Noco > 1) error("not prepared for Noco= %d", Noco);
        } // iAtom

        std::vector<double> Veff(nblocks*n4x4x4, 0.);
        double constexpr scale_V = 1.0;
        set(Veff.data(), nblocks*size_t(64), Vtot, scale_V);
        stat += green_function::update_potential(plan, pg.grid_blocks(), Veff, AtomMatrices, echo, Noco);

        view2D<double> kpoint_mesh;
        // get a kpoint mesh controlled by +hamiltonian.kmesh.x .y .z, the same for each energy point
        auto const nkpoints = brillouin_zone::get_kpoint_mesh(kpoint_mesh);

        std::vector<Complex> energy_weights;
        auto const energy_mesh = get_energy_mesh(energy_weights, echo);
        int const nEpoints = energy_mesh.size();
        if (echo > 0) std::printf("# energy_contour::integration with %d energy points and %d k-points\n", nEpoints, nkpoints);

        Complex constexpr zero = 0;
        view2D<Complex> rho_c(nblocks, n4x4x4, zero); // complex density

        for (int iEpoint{0}; iEpoint < nEpoints; ++iEpoint) {
            Complex const energy_weight = energy_weights[iEpoint];

            Complex const energy = energy_mesh.at(iEpoint) + Fermi_level;
            char energy_parameter_label[64];
            std::snprintf(energy_parameter_label, 64, "(%g %s, %g %s)", (energy.real() - Fermi_level)*eV, _eV, energy.imag()*Kelvin, _Kelvin);
            if (echo > 5) std::printf("# energy parameter %s with weight (%g, %g)\n", energy_parameter_label, energy_weight.real(), energy_weight.imag());

            stat += green_function::update_energy_parameter(plan, energy, dV, echo, Noco);

            view2D<Complex> rho_E(nblocks, n4x4x4, zero);

            for (int ikpoint{0}; ikpoint < nkpoints; ++ikpoint) {
                double const *const kpoint = kpoint_mesh[ikpoint];
                Complex const kpoint_weight = kpoint[brillouin_zone::WEIGHT];

                if (echo + check > 5) std::printf("# solve Green function for E=%s, k-point=[%g %g %g]\n",
                                                energy_parameter_label, kpoint[0], kpoint[1], kpoint[2]);
                if (0 == check) {
                    stat += green_function::update_phases(plan, kpoint, echo, Noco);

                    view2D<Complex> rho_Ek(nblocks, n4x4x4, zero);

                    stat += solver_->solve(rho_Ek[0], nblocks, iterations, echo);

                    add_product(rho_E[0], nblocks*n4x4x4, rho_Ek[0], kpoint_weight); // accumulate density over k-points
                    auto const rho_integral = mpi_parallel::sum(sum(rho_Ek[0], nblocks*n4x4x4).imag(), comm)*dV;
                    if (echo + check > 6) std::printf("# solved Green function for E=%s, k-point=[%g %g %g] has %g electrons\n",
                                                    energy_parameter_label, kpoint[0], kpoint[1], kpoint[2], rho_integral);
                } // check

            } // ikpoint
            if (0 == check) {
                add_product(rho_c[0], nblocks*n4x4x4, rho_E[0], energy_weight); // accumulate density over E-points
                auto const rho_integral = mpi_parallel::sum(sum(rho_E[0], nblocks*n4x4x4).imag(), comm)*dV;
                if (echo + check > 3) std::printf("# solved Green function for E=%s has %g electrons\n",
                                                                 energy_parameter_label, rho_integral);
            } // check
            if (echo + check > 2) std::printf("# solved Green function for E=%s\n", energy_parameter_label);
        } // iEpoint


        view2D<double> rho(nblocks, n4x4x4, 0.0);
        for (uint32_t ib{0}; ib < nblocks; ++ib) {
            for (int i444{0}; i444 < 64; ++i444) {
                rho(ib,i444) = rho_c(ib,i444).real();
            } // i444
        } // ib

        {
            auto const rho_integral = mpi_parallel::sum(sum(rho[0], nblocks*n4x4x4), comm)*dV;
            if (echo + check > 3) std::printf("# solved density has %g electrons\n", rho_integral);
            if (echo > 3) std::printf("# rank#%i maxval rho= %g a.u.\n", me, maxval(rho[0], nblocks*n4x4x4));
        }

        // interpolation density from 4*4*4 to 8*8*8 block could be done here
        if (echo > 3) std::printf("# interpolate density from 4x4x4 to 8x8x8\n");
        parallel_poisson::block_interpolation(rho_888, rho[0], pg, echo, 1., "density");

        {
            auto const rho_integral = mpi_parallel::sum(sum(rho_888, nblocks*size_t(8*8*8)), comm)*dV;
            if (echo + check > 3) std::printf("# interpolated density has %g electrons\n", rho_integral);
            if (echo > 3) std::printf("# rank#%i maxval rho= %g a.u.\n", me, maxval(rho_888, nblocks*512));
        }

        if (echo > 3) std::printf("# density integrated over %d energy points\n", nEpoints);
        return stat;
    } // integrate




    int Gauss_Legendre_quadrature(double x[], double w[], unsigned const number, int const echo) {
        // from juKKR: https://iffgit.fz-juelich.de/kkr/jukkr/-/blob/master/source/KKRhost/gauleg.f90
        int n = std::max(number, 1u);
        if (n > 64) { n = ((n - 1)/8 + 1)*8; } else
        if (n > 32) { n = ((n - 1)/4 + 1)*4; }
        if (number != n && echo > 5) std::printf("# %s(%d --> %d)\n", __func__, number, n);

#define xi(INDEX) x[(INDEX) - 1]
#define wi(INDEX) w[(INDEX) - 1]

#if 1
    if (number==0) {
        return 0; // ok
    } else
    if (n==1) {
      xi(1) = 0.e0;
      wi(1) = 2.e0;
    } else
    if (n==2) {
      xi(1) = -57735026918962576451.e-20;
      wi(1) = 10000000000000000000.e-19;
    } else
    if (n==3) {
      xi(1) = -77459666924148337704.e-20;
      xi(2) = 0.e0;
      wi(1) = 55555555555555555556.e-20;
      wi(2) = 88888888888888888889.e-20;
    } else 
    if (n==4) {
      xi(1) = -86113631159405257522.e-20;
      xi(2) = -33998104358485626480.e-20;
      wi(1) = 34785484513745385737.e-20;
      wi(2) = 65214515486254614263.e-20;
      
    }

    else if (n==5) {
      xi(1) = -90617984593866399280.e-20;
      xi(2) = -53846931010568309104.e-20;
      xi(3) = 0.e0;
      wi(1) = 23692688505618908751.e-20;
      wi(2) = 47862867049936646804.e-20;
      wi(3) = 56888888888888888889.e-20;
      
    }

    else if (n==6) {
      xi(1) = -93246951420315202781.e-20;
      xi(2) = -66120938646626451366.e-20;
      xi(3) = -23861918608319690863.e-20;
      wi(1) = 17132449237917034504.e-20;
      wi(2) = 36076157304813860757.e-20;
      wi(3) = 46791393457269104739.e-20;
      
    }

    else if (n==7) {
      xi(1) = -94910791234275852453.e-20;
      xi(2) = -74153118559939443986.e-20;
      xi(3) = -40584515137739716691.e-20;
      xi(4) = 0.e0;
      wi(1) = 12948496616886969327.e-20;
      wi(2) = 27970539148927666790.e-20;
      wi(3) = 38183005050511894495.e-20;
      wi(4) = 41795918367346938776.e-20;
      
    }

    else if (n==8) {
      xi(1) = -96028985649753623168.e-20;
      xi(2) = -79666647741362673959.e-20;
      xi(3) = -52553240991632898582.e-20;
      xi(4) = -18343464249564980494.e-20;
      wi(1) = 10122853629037625915.e-20;
      wi(2) = 22238103445337447054.e-20;
      wi(3) = 31370664587788728734.e-20;
      wi(4) = 36268378337836198297.e-20;
      

    }

    else if (n==9) {
      xi(1) = -96816023950762608984.e-20;
      xi(2) = -83603110732663579430.e-20;
      xi(3) = -61337143270059039731.e-20;
      xi(4) = -32425342340380892904.e-20;
      xi(5) = 0.e0;
      wi(1) = 81274388361574411972.e-21;
      wi(2) = 18064816069485740406.e-20;
      wi(3) = 26061069640293546232.e-20;
      wi(4) = 31234707704000284007.e-20;
      wi(5) = 33023935500125976316.e-20;
      
    }

    else if (n==10) {
      xi(1) = -97390652851717172008.e-20;
      xi(2) = -86506336668898451073.e-20;
      xi(3) = -67940956829902440623.e-20;
      xi(4) = -43339539412924719080.e-20;
      xi(5) = -14887433898163121088.e-20;
      wi(1) = 66671344308688137594.e-21;
      wi(2) = 14945134915058059315.e-20;
      wi(3) = 21908636251598204400.e-20;
      wi(4) = 26926671930999635509.e-20;
      wi(5) = 29552422471475287017.e-20;
      
    }

    else if (n==11) {
      xi(1) = -97822865814605699280.e-20;
      xi(2) = -88706259976809529908.e-20;
      xi(3) = -73015200557404932409.e-20;
      xi(4) = -51909612920681181593.e-20;
      xi(5) = -26954315595234497233.e-20;
      xi(6) = 0.e0;
      wi(1) = 55668567116173666483.e-21;
      wi(2) = 12558036946490462463.e-20;
      wi(3) = 18629021092773425143.e-20;
      wi(4) = 23319376459199047992.e-20;
      wi(5) = 26280454451024666218.e-20;
      wi(6) = 27292508677790063071.e-20;
      
    }

    else if (n==12) {
      xi(1) = -98156063424671925069.e-20;
      xi(2) = -90411725637047485668.e-20;
      xi(3) = -76990267419430468704.e-20;
      xi(4) = -58731795428661744730.e-20;
      xi(5) = -36783149899818019375.e-20;
      xi(6) = -12523340851146891547.e-20;
      wi(1) = 47175336386511827195.e-21;
      wi(2) = 10693932599531843096.e-20;
      wi(3) = 16007832854334622633.e-20;
      wi(4) = 20316742672306592175.e-20;
      wi(5) = 23349253653835480876.e-20;
      wi(6) = 24914704581340278500.e-20;
      
    }

    else if (n==13) {
      xi(1) = -98418305471858814947.e-20;
      xi(2) = -91759839922297796521.e-20;
      xi(3) = -80157809073330991279.e-20;
      xi(4) = -64234933944034022064.e-20;
      xi(5) = -44849275103644685288.e-20;
      xi(6) = -23045831595513479407.e-20;
      xi(7) = 0.e0;
      wi(1) = 40484004765315879520.e-21;
      wi(2) = 92121499837728447914.e-21;
      wi(3) = 13887351021978723846.e-20;
      wi(4) = 17814598076194573828.e-20;
      wi(5) = 20781604753688850231.e-20;
      wi(6) = 22628318026289723841.e-20;
      wi(7) = 23255155323087391019.e-20;
      
    }

    else if (n==14) {
      xi(1) = -98628380869681233884.e-20;
      xi(2) = -92843488366357351734.e-20;
      xi(3) = -82720131506976499319.e-20;
      xi(4) = -68729290481168547015.e-20;
      xi(5) = -51524863635815409197.e-20;
      xi(6) = -31911236892788976044.e-20;
      xi(7) = -10805494870734366207.e-20;
      wi(1) = 35119460331751863032.e-21;
      wi(2) = 80158087159760209806.e-21;
      wi(3) = 12151857068790318469.e-20;
      wi(4) = 15720316715819353457.e-20;
      wi(5) = 18553839747793781374.e-20;
      wi(6) = 20519846372129560397.e-20;
      wi(7) = 21526385346315779020.e-20;
      
    }

    else if (n==15) {
      xi(1) = -98799251802048542849.e-20;
      xi(2) = -93727339240070590431.e-20;
      xi(3) = -84820658341042721620.e-20;
      xi(4) = -72441773136017004742.e-20;
      xi(5) = -57097217260853884754.e-20;
      xi(6) = -39415134707756336990.e-20;
      xi(7) = -20119409399743452230.e-20;
      xi(8) = 0.e0;
      wi(1) = 30753241996117268355.e-21;
      wi(2) = 70366047488108124709.e-21;
      wi(3) = 10715922046717193501.e-20;
      wi(4) = 13957067792615431445.e-20;
      wi(5) = 16626920581699393355.e-20;
      wi(6) = 18616100001556221103.e-20;
      wi(7) = 19843148532711157646.e-20;
      wi(8) = 20257824192556127288.e-20;
      
    }

    else if (n==16) {
      xi(1) = -98940093499164993260.e-20;
      xi(2) = -94457502307323257608.e-20;
      xi(3) = -86563120238783174388.e-20;
      xi(4) = -75540440835500303390.e-20;
      xi(5) = -61787624440264374845.e-20;
      xi(6) = -45801677765722738634.e-20;
      xi(7) = -28160355077925891323.e-20;
      xi(8) = -95012509837637440185.e-21;
      wi(1) = 27152459411754094852.e-21;
      wi(2) = 62253523938647892863.e-21;
      wi(3) = 95158511682492784810.e-21;
      wi(4) = 12462897125553387205.e-20;
      wi(5) = 14959598881657673208.e-20;
      wi(6) = 16915651939500253819.e-20;
      wi(7) = 18260341504492358887.e-20;
      wi(8) = 18945061045506849629.e-20;
      
    }

    else if (n==17) {
      xi(1) = -99057547531441733568.e-20;
      xi(2) = -95067552176876776122.e-20;
      xi(3) = -88023915372698590212.e-20;
      xi(4) = -78151400389680140693.e-20;
      xi(5) = -65767115921669076585.e-20;
      xi(6) = -51269053708647696789.e-20;
      xi(7) = -35123176345387631530.e-20;
      xi(8) = -17848418149584785585.e-20;
      xi(9) = 0.e0;
      wi(1) = 24148302868547931960.e-21;
      wi(2) = 55459529373987201129.e-21;
      wi(3) = 85036148317179180884.e-21;
      wi(4) = 11188384719340397109.e-20;
      wi(5) = 13513636846852547329.e-20;
      wi(6) = 15404576107681028808.e-20;
      wi(7) = 16800410215645004451.e-20;
      wi(8) = 17656270536699264633.e-20;
      wi(9) = 17944647035620652546.e-20;
      
    }

    else if (n==18) {
      xi(1) = -99156516842093094673.e-20;
      xi(2) = -95582394957139775518.e-20;
      xi(3) = -89260246649755573921.e-20;
      xi(4) = -80370495897252311568.e-20;
      xi(5) = -69168704306035320787.e-20;
      xi(6) = -55977083107394753461.e-20;
      xi(7) = -41175116146284264604.e-20;
      xi(8) = -25188622569150550959.e-20;
      xi(9) = -84775013041735301242.e-21;
      wi(1) = 21616013526483310313.e-21;
      wi(2) = 49714548894969796453.e-21;
      wi(3) = 76425730254889056529.e-21;
      wi(4) = 10094204410628716556.e-20;
      wi(5) = 12255520671147846018.e-20;
      wi(6) = 14064291467065065120.e-20;
      wi(7) = 15468467512626524493.e-20;
      wi(8) = 16427648374583272299.e-20;
      wi(9) = 16914238296314359184.e-20;
      
    }

    else if (n==19) {
      xi(1) = -99240684384358440319.e-20;
      xi(2) = -96020815213483003085.e-20;
      xi(3) = -90315590361481790164.e-20;
      xi(4) = -82271465653714282498.e-20;
      xi(5) = -72096617733522937862.e-20;
      xi(6) = -60054530466168102347.e-20;
      xi(7) = -46457074137596094572.e-20;
      xi(8) = -31656409996362983199.e-20;
      xi(9) = -16035864564022537587.e-20;
      xi(10) = 0.e0;
      wi(1) = 19461788229726477036.e-21;
      wi(2) = 44814226765699600333.e-21;
      wi(3) = 69044542737641226581.e-21;
      wi(4) = 91490021622449999464.e-21;
      wi(5) = 11156664554733399472.e-20;
      wi(6) = 12875396253933622768.e-20;
      wi(7) = 14260670217360661178.e-20;
      wi(8) = 15276604206585966678.e-20;
      wi(9) = 15896884339395434765.e-20;
      wi(10) = 16105444984878369598.e-20;
      
    }

    else if (n==20) {
      xi(1) = -99312859918509492479.e-20;
      xi(2) = -96397192727791379127.e-20;
      xi(3) = -91223442825132590587.e-20;
      xi(4) = -83911697182221882339.e-20;
      xi(5) = -74633190646015079261.e-20;
      xi(6) = -63605368072651502545.e-20;
      xi(7) = -51086700195082709800.e-20;
      xi(8) = -37370608871541956067.e-20;
      xi(9) = -22778585114164507808.e-20;
      xi(10) = -76526521133497333755.e-21;
      wi(1) = 17614007139152118312.e-21;
      wi(2) = 40601429800386941331.e-21;
      wi(3) = 62672048334109063570.e-21;
      wi(4) = 83276741576704748725.e-21;
      wi(5) = 10193011981724043504.e-20;
      wi(6) = 11819453196151841731.e-20;
      wi(7) = 13168863844917662690.e-20;
      wi(8) = 14209610931838205133.e-20;
      wi(9) = 14917298647260374679.e-20;
      wi(10) = 15275338713072585070.e-20;
      
    }

    else if (n==21) {
      xi(1) = -99375217062038950026.e-20;
      xi(2) = -96722683856630629432.e-20;
      xi(3) = -92009933415040082879.e-20;
      xi(4) = -85336336458331728365.e-20;
      xi(5) = -76843996347567790862.e-20;
      xi(6) = -66713880419741231931.e-20;
      xi(7) = -55161883588721980706.e-20;
      xi(8) = -42434212020743878357.e-20;
      xi(9) = -28802131680240109660.e-20;
      xi(10) = -14556185416089509094.e-20;
      xi(11) = 0.e0;
      wi(1) = 16017228257774333324.e-21;
      wi(2) = 36953789770852493800.e-21;
      wi(3) = 57134425426857208284.e-21;
      wi(4) = 76100113628379302017.e-21;
      wi(5) = 93444423456033861553.e-21;
      wi(6) = 10879729916714837766.e-20;
      wi(7) = 12183141605372853420.e-20;
      wi(8) = 13226893863333746178.e-20;
      wi(9) = 13988739479107315472.e-20;
      wi(10) = 14452440398997005906.e-20;
      wi(11) = 14608113364969042719.e-20;
      
    }

    else if (n==22) {
      xi(1) = -99429458548239929207.e-20;
      xi(2) = -97006049783542872712.e-20;
      xi(3) = -92695677218717400052.e-20;
      xi(4) = -86581257772030013654.e-20;
      xi(5) = -78781680597920816200.e-20;
      xi(6) = -69448726318668278005.e-20;
      xi(7) = -58764040350691159296.e-20;
      xi(8) = -46935583798675702641.e-20;
      xi(9) = -34193582089208422516.e-20;
      xi(10) = -20786042668822128548.e-20;
      xi(11) = -69739273319722221214.e-21;
      wi(1) = 14627995298272200685.e-21;
      wi(2) = 33774901584814154793.e-21;
      wi(3) = 52293335152683285940.e-21;
      wi(4) = 69796468424520488095.e-21;
      wi(5) = 85941606217067727414.e-21;
      wi(6) = 10041414444288096493.e-20;
      wi(7) = 11293229608053921839.e-20;
      wi(8) = 12325237681051242429.e-20;
      wi(9) = 13117350478706237073.e-20;
      wi(10) = 13654149834601517135.e-20;
      wi(11) = 13925187285563199338.e-20;
      
    }

    else if (n==23) {
      xi(1) = -99476933499755212352.e-20;
      xi(2) = -97254247121811523196.e-20;
      xi(3) = -93297108682601610235.e-20;
      xi(4) = -87675235827044166738.e-20;
      xi(5) = -80488840161883989215.e-20;
      xi(6) = -71866136313195019446.e-20;
      xi(7) = -61960987576364615639.e-20;
      xi(8) = -50950147784600754969.e-20;
      xi(9) = -39030103803029083142.e-20;
      xi(10) = -26413568097034493053.e-20;
      xi(11) = -13325682429846611093.e-20;
      xi(12) = 0.e0;
      wi(1) = 13411859487141772081.e-21;
      wi(2) = 30988005856979444311.e-21;
      wi(3) = 48037671731084668572.e-21;
      wi(4) = 64232421408525852127.e-21;
      wi(5) = 79281411776718954923.e-21;
      wi(6) = 92915766060035147477.e-21;
      wi(7) = 10489209146454141007.e-20;
      wi(8) = 11499664022241136494.e-20;
      wi(9) = 12304908430672953047.e-20;
      wi(10) = 12890572218808214998.e-20;
      wi(11) = 13246203940469661737.e-20;
      wi(12) = 13365457218610617535.e-20;
      
    }

    else if (n==24) {
      xi(1) = -99518721999702136018.e-20;
      xi(2) = -97472855597130949820.e-20;
      xi(3) = -93827455200273275852.e-20;
      xi(4) = -88641552700440103421.e-20;
      xi(5) = -82000198597390292195.e-20;
      xi(6) = -74012419157855436424.e-20;
      xi(7) = -64809365193697556925.e-20;
      xi(8) = -54542147138883953566.e-20;
      xi(9) = -43379350762604513849.e-20;
      xi(10) = -31504267969616337439.e-20;
      xi(11) = -19111886747361630916.e-20;
      xi(12) = -64056892862605626085.e-21;
      wi(1) = 12341229799987199547.e-21;
      wi(2) = 28531388628933663181.e-21;
      wi(3) = 44277438817419806169.e-21;
      wi(4) = 59298584915436780746.e-21;
      wi(5) = 73346481411080305734.e-21;
      wi(6) = 86190161531953275917.e-21;
      wi(7) = 97618652104113888270.e-21;
      wi(8) = 10744427011596563478.e-20;
      wi(9) = 11550566805372560135.e-20;
      wi(10) = 12167047292780339120.e-20;
      wi(11) = 12583745634682829612.e-20;
      wi(12) = 12793819534675215697.e-20;
      
    }

    else if (n==25) {
      xi(1) = -99555696979049809791.e-20;
      xi(2) = -97666392145951751150.e-20;
      xi(3) = -94297457122897433941.e-20;
      xi(4) = -89499199787827536885.e-20;
      xi(5) = -83344262876083400142.e-20;
      xi(6) = -75925926303735763058.e-20;
      xi(7) = -67356636847346836449.e-20;
      xi(8) = -57766293024122296772.e-20;
      xi(9) = -47300273144571496052.e-20;
      xi(10) = -36117230580938783774.e-20;
      xi(11) = -24386688372098843205.e-20;
      xi(12) = -12286469261071039639.e-20;
      xi(13) = 0.e0;
      wi(1) = 11393798501026287948.e-21;
      wi(2) = 26354986615032137262.e-21;
      wi(3) = 40939156701306312656.e-21;
      wi(4) = 54904695975835191926.e-21;
      wi(5) = 68038333812356917207.e-21;
      wi(6) = 80140700335001018013.e-21;
      wi(7) = 91028261982963649811.e-21;
      wi(8) = 10053594906705064420.e-20;
      wi(9) = 10851962447426365312.e-20;
      wi(10) = 11485825914571164834.e-20;
      wi(11) = 11945576353578477223.e-20;
      wi(12) = 12224244299031004169.e-20;
      wi(13) = 12317605372671545120.e-20;
      
    }

    else if (n==26) {
      xi(1) = -99588570114561692900.e-20;
      xi(2) = -97838544595647099110.e-20;
      xi(3) = -94715906666171425014.e-20;
      xi(4) = -90263786198430707422.e-20;
      xi(5) = -84544594278849801880.e-20;
      xi(6) = -77638594882067885619.e-20;
      xi(7) = -69642726041995726486.e-20;
      xi(8) = -60669229301761806323.e-20;
      xi(9) = -50844071482450571770.e-20;
      xi(10) = -40305175512348630648.e-20;
      xi(11) = -29200483948595689514.e-20;
      xi(12) = -17685882035689018397.e-20;
      xi(13) = -59230093429313207094.e-21;
      wi(1) = 10551372617343007156.e-21;
      wi(2) = 24417851092631908790.e-21;
      wi(3) = 37962383294362763950.e-21;
      wi(4) = 50975825297147811998.e-21;
      wi(5) = 63274046329574835539.e-21;
      wi(6) = 74684149765659745887.e-21;
      wi(7) = 85045894313485239210.e-21;
      wi(8) = 94213800355914148464.e-21;
      wi(9) = 10205916109442542324.e-20;
      wi(10) = 10847184052857659066.e-20;
      wi(11) = 11336181654631966655.e-20;
      wi(12) = 11666044348529658204.e-20;
      wi(13) = 11832141527926227652.e-20;
      
    }

    else if (n==27) {
      xi(1) = -99617926288898856694.e-20;
      xi(2) = -97992347596150122286.e-20;
      xi(3) = -95090055781470500685.e-20;
      xi(4) = -90948232067749110430.e-20;
      xi(5) = -85620790801829449030.e-20;
      xi(6) = -79177163907050822714.e-20;
      xi(7) = -71701347373942369929.e-20;
      xi(8) = -63290797194649514093.e-20;
      xi(9) = -54055156457945689490.e-20;
      xi(10) = -44114825175002688059.e-20;
      xi(11) = -33599390363850889973.e-20;
      xi(12) = -22645936543953685886.e-20;
      xi(13) = -11397258560952996693.e-20;
      xi(14) = 0.e0;
      wi(1) = 97989960512943602612.e-22;
      wi(2) = 22686231596180623196.e-21;
      wi(3) = 35297053757419711023.e-21;
      wi(4) = 47449412520615062704.e-21;
      wi(5) = 58983536859833599110.e-21;
      wi(6) = 69748823766245592984.e-21;
      wi(7) = 79604867773057771263.e-21;
      wi(8) = 88423158543756950194.e-21;
      wi(9) = 96088727370028507566.e-21;
      wi(10) = 10250163781774579867.e-20;
      wi(11) = 10757828578853318721.e-20;
      wi(12) = 11125248835684519267.e-20;
      wi(13) = 11347634610896514862.e-20;
      wi(14) = 11422086737895698905.e-20;
      
    }

    else if (n==28) {
      xi(1) = -99644249757395444995.e-20;
      xi(2) = -98130316537087275369.e-20;
      xi(3) = -95425928062893819725.e-20;
      xi(4) = -91563302639213207387.e-20;
      xi(5) = -86589252257439504894.e-20;
      xi(6) = -80564137091717917145.e-20;
      xi(7) = -73561087801363177203.e-20;
      xi(8) = -65665109403886496122.e-20;
      xi(9) = -56972047181140171931.e-20;
      xi(10) = -47587422495511826103.e-20;
      xi(11) = -37625151608907871022.e-20;
      xi(12) = -27206162763517807768.e-20;
      xi(13) = -16456928213338077128.e-20;
      xi(14) = -55079289884034270427.e-21;
      wi(1) = 91242825930945177388.e-22;
      wi(2) = 21132112592771259752.e-21;
      wi(3) = 32901427782304379978.e-21;
      wi(4) = 44272934759004227840.e-21;
      wi(5) = 55107345675716745431.e-21;
      wi(6) = 65272923966999595793.e-21;
      wi(7) = 74646214234568779024.e-21;
      wi(8) = 83113417228901218390.e-21;
      wi(9) = 90571744393032840942.e-21;
      wi(10) = 96930657997929915850.e-21;
      wi(11) = 10211296757806076981.e-20;
      wi(12) = 10605576592284641791.e-20;
      wi(13) = 10871119225829413525.e-20;
      wi(14) = 11004701301647519628.e-20;
      
    }

    else if (n==29) {
      xi(1) = -99667944226059658616.e-20;
      xi(2) = -98254550526141317487.e-20;
      xi(3) = -95728559577808772580.e-20;
      xi(4) = -92118023295305878509.e-20;
      xi(5) = -87463780492010279042.e-20;
      xi(6) = -81818548761525244499.e-20;
      xi(7) = -75246285173447713391.e-20;
      xi(8) = -67821453760268651516.e-20;
      xi(9) = -59628179713822782038.e-20;
      xi(10) = -50759295512422764210.e-20;
      xi(11) = -41315288817400866389.e-20;
      xi(12) = -31403163786763993495.e-20;
      xi(13) = -21135228616600107451.e-20;
      xi(14) = -10627823013267923017.e-20;
      xi(15) = 0.e0;
      wi(1) = 85169038787464096543.e-22;
      wi(2) = 19732085056122705984.e-21;
      wi(3) = 30740492202093622644.e-21;
      wi(4) = 41402062518682836105.e-21;
      wi(5) = 51594826902497923913.e-21;
      wi(6) = 61203090657079138542.e-21;
      wi(7) = 70117933255051278570.e-21;
      wi(8) = 78238327135763783828.e-21;
      wi(9) = 85472257366172527545.e-21;
      wi(10) = 91737757139258763348.e-21;
      wi(11) = 96963834094408606302.e-21;
      wi(12) = 10109127375991496612.e-20;
      wi(13) = 10407331007772937391.e-20;
      wi(14) = 10587615509732094141.e-20;
      wi(15) = 10647938171831424425.e-20;
      
    }

    else if (n==30) {
      xi(1) = -99689348407464954027.e-20;
      xi(2) = -98366812327974720997.e-20;
      xi(3) = -96002186496830751222.e-20;
      xi(4) = -92620004742927432588.e-20;
      xi(5) = -88256053579205268154.e-20;
      xi(6) = -82956576238276839744.e-20;
      xi(7) = -76777743210482619492.e-20;
      xi(8) = -69785049479331579693.e-20;
      xi(9) = -62052618298924286114.e-20;
      xi(10) = -53662414814201989926.e-20;
      xi(11) = -44703376953808917678.e-20;
      xi(12) = -35270472553087811347.e-20;
      xi(13) = -25463692616788984644.e-20;
      xi(14) = -15386991360858354696.e-20;
      xi(15) = -51471842555317695833.e-21;
      wi(1) = 79681924961666056155.e-22;
      wi(2) = 18466468311090959142.e-21;
      wi(3) = 28784707883323369350.e-21;
      wi(4) = 38799192569627049597.e-21;
      wi(5) = 48402672830594052903.e-21;
      wi(6) = 57493156217619066482.e-21;
      wi(7) = 65974229882180495128.e-21;
      wi(8) = 73755974737705206268.e-21;
      wi(9) = 80755895229420215355.e-21;
      wi(10) = 86899787201082979802.e-21;
      wi(11) = 92122522237786128718.e-21;
      wi(12) = 96368737174644259639.e-21;
      wi(13) = 99593420586795267063.e-21;
      wi(14) = 10176238974840550460.e-20;
      wi(15) = 10285265289355884034.e-20;
      
    }

    else if (n==31) {
      xi(1) = -99708748181947707406.e-20;
      xi(2) = -98468590966515248400.e-20;
      xi(3) = -96250392509294966179.e-20;
      xi(4) = -93075699789664816496.e-20;
      xi(5) = -88976002994827104337.e-20;
      xi(6) = -83992032014626734009.e-20;
      xi(7) = -78173314841662494041.e-20;
      xi(8) = -71577678458685328391.e-20;
      xi(9) = -64270672292426034618.e-20;
      xi(10) = -56324916140714926272.e-20;
      xi(11) = -47819378204490248044.e-20;
      xi(12) = -38838590160823294306.e-20;
      xi(13) = -29471806998170161662.e-20;
      xi(14) = -19812119933557062877.e-20;
      xi(15) = -99555312152341520325.e-21;
      xi(16) = 0.e0;
      wi(1) = 74708315792487758587.e-22;
      wi(2) = 17318620790310582463.e-21;
      wi(3) = 27009019184979421801.e-21;
      wi(4) = 36432273912385464024.e-21;
      wi(5) = 45493707527201102902.e-21;
      wi(6) = 54103082424916853712.e-21;
      wi(7) = 62174786561028426910.e-21;
      wi(8) = 69628583235410366168.e-21;
      wi(9) = 76390386598776616426.e-21;
      wi(10) = 82392991761589263904.e-21;
      wi(11) = 87576740608477876126.e-21;
      wi(12) = 91890113893641478215.e-21;
      wi(13) = 95290242912319512807.e-21;
      wi(14) = 97743335386328725093.e-21;
      wi(15) = 99225011226672307875.e-21;
      wi(16) = 99720544793426451428.e-21;
      
    }

    else if (n==32) {
      xi(1) = -99726386184948156354.e-20;
      xi(2) = -98561151154526833540.e-20;
      xi(3) = -96476225558750643077.e-20;
      xi(4) = -93490607593773968917.e-20;
      xi(5) = -89632115576605212397.e-20;
      xi(6) = -84936761373256997013.e-20;
      xi(7) = -79448379596794240696.e-20;
      xi(8) = -73218211874028968039.e-20;
      xi(9) = -66304426693021520098.e-20;
      xi(10) = -58771575724076232904.e-20;
      xi(11) = -50689990893222939002.e-20;
      xi(12) = -42135127613063534536.e-20;
      xi(13) = -33186860228212764978.e-20;
      xi(14) = -23928736225213707454.e-20;
      xi(15) = -14447196158279649349.e-20;
      xi(16) = -48307665687738316235.e-21;
      wi(1) = 70186100094700966004.e-22;
      wi(2) = 16274394730905670605.e-21;
      wi(3) = 25392065309262059456.e-21;
      wi(4) = 34273862913021433103.e-21;
      wi(5) = 42835898022226680657.e-21;
      wi(6) = 50998059262376176196.e-21;
      wi(7) = 58684093478535547145.e-21;
      wi(8) = 65822222776361846838.e-21;
      wi(9) = 72345794108848506225.e-21;
      wi(10) = 78193895787070306472.e-21;
      wi(11) = 83311924226946755222.e-21;
      wi(12) = 87652093004403811143.e-21;
      wi(13) = 91173878695763884713.e-21;
      wi(14) = 93844399080804565639.e-21;
      wi(15) = 95638720079274859419.e-21;
      wi(16) = 96540088514727800567.e-21;
      
    }

    else if (n==36) {
      xi(1) = -99783046248408583620.e-20;
      xi(2) = -98858647890221223807.e-20;
      xi(3) = -97202769104969794934.e-20;
      xi(4) = -94827298439950754520.e-20;
      xi(5) = -91749777451565906608.e-20;
      xi(6) = -87992980089039713198.e-20;
      xi(7) = -83584716699247530642.e-20;
      xi(8) = -78557623013220651283.e-20;
      xi(9) = -72948917159355658209.e-20;
      xi(10) = -66800123658552106210.e-20;
      xi(11) = -60156765813598053508.e-20;
      xi(12) = -53068028592624516164.e-20;
      xi(13) = -45586394443342026721.e-20;
      xi(14) = -37767254711968921632.e-20;
      xi(15) = -29668499534402827050.e-20;
      xi(16) = -21350089231686557894.e-20;
      xi(17) = -12873610380938478865.e-20;
      xi(18) = -43018198473708607227.e-21;
      wi(1) = 55657196642450453613.e-22;
      wi(2) = 12915947284065574405.e-21;
      wi(3) = 20181515297735471532.e-21;
      wi(4) = 27298621498568779094.e-21;
      wi(5) = 34213810770307229921.e-21;
      wi(6) = 40875750923644895474.e-21;
      wi(7) = 47235083490265978417.e-21;
      wi(8) = 53244713977759919092.e-21;
      wi(9) = 58860144245324817310.e-21;
      wi(10) = 64039797355015489556.e-21;
      wi(11) = 68745323835736442614.e-21;
      wi(12) = 72941885005653061354.e-21;
      wi(13) = 76598410645870674529.e-21;
      wi(14) = 79687828912071601909.e-21;
      wi(15) = 82187266704339709517.e-21;
      wi(16) = 84078218979661934933.e-21;
      wi(17) = 85346685739338627492.e-21;
      wi(18) = 85983275670394747490.e-21;
      
    }

    else if (n==40) {
      xi(1) = -99823770971055920035.e-20;
      xi(2) = -99072623869945700645.e-20;
      xi(3) = -97725994998377426266.e-20;
      xi(4) = -95791681921379165580.e-20;
      xi(5) = -93281280827867653336.e-20;
      xi(6) = -90209880696887429673.e-20;
      xi(7) = -86595950321225950382.e-20;
      xi(8) = -82461223083331166320.e-20;
      xi(9) = -77830565142651938769.e-20;
      xi(10) = -72731825518992710328.e-20;
      xi(11) = -67195668461417954838.e-20;
      xi(12) = -61255388966798023795.e-20;
      xi(13) = -54946712509512820208.e-20;
      xi(14) = -48307580168617871291.e-20;
      xi(15) = -41377920437160500152.e-20;
      xi(16) = -34199409082575847301.e-20;
      xi(17) = -26815218500725368114.e-20;
      xi(18) = -19269758070137109972.e-20;
      xi(19) = -11608407067525520848.e-20;
      xi(20) = -38772417506050821933.e-21;
      wi(1) = 45212770985331912585.e-22;
      wi(2) = 10498284531152813615.e-21;
      wi(3) = 16421058381907888713.e-21;
      wi(4) = 22245849194166957262.e-21;
      wi(5) = 27937006980023401098.e-21;
      wi(6) = 33460195282547847393.e-21;
      wi(7) = 38782167974472017640.e-21;
      wi(8) = 43870908185673271992.e-21;
      wi(9) = 48695807635072232061.e-21;
      wi(10) = 53227846983936824355.e-21;
      wi(11) = 57439769099391551367.e-21;
      wi(12) = 61306242492928939167.e-21;
      wi(13) = 64804013456601038075.e-21;
      wi(14) = 67912045815233903826.e-21;
      wi(15) = 70611647391286779695.e-21;
      wi(16) = 72886582395804059061.e-21;
      wi(17) = 74723169057968264200.e-21;
      wi(18) = 76110361900626242372.e-21;
      wi(19) = 77039818164247965588.e-21;
      wi(20) = 77505947978424811264.e-21;
      
    }

    else if (n==44) {
      xi(1) = -99854020063677422494.e-20;
      xi(2) = -99231639213851580848.e-20;
      xi(3) = -98115183307791396666.e-20;
      xi(4) = -96509965042249313939.e-20;
      xi(5) = -94423950911819409920.e-20;
      xi(6) = -91867525998417577432.e-20;
      xi(7) = -88853423828604320234.e-20;
      xi(8) = -85396659500471037873.e-20;
      xi(9) = -81514453964513501049.e-20;
      xi(10) = -77226147924875589902.e-20;
      xi(11) = -72553105366071700261.e-20;
      xi(12) = -67518607066612236533.e-20;
      xi(13) = -62147734590357584780.e-20;
      xi(14) = -56467245318547076842.e-20;
      xi(15) = -50505439138820231798.e-20;
      xi(16) = -44292017452541148383.e-20;
      xi(17) = -37857935201470713251.e-20;
      xi(18) = -31235246650278581224.e-20;
      xi(19) = -24456945692820125151.e-20;
      xi(20) = -17556801477551678575.e-20;
      xi(21) = -10569190170865324712.e-20;
      xi(22) = -35289236964135359058.e-21;
      wi(1) = 37454048031127775152.e-22;
      wi(2) = 87004813675248441226.e-22;
      wi(3) = 13619586755579985520.e-21;
      wi(4) = 18471481736814749172.e-21;
      wi(5) = 23231481902019210629.e-21;
      wi(6) = 27875782821281010081.e-21;
      wi(7) = 32381222812069820881.e-21;
      wi(8) = 36725347813808873643.e-21;
      wi(9) = 40886512310346218908.e-21;
      wi(10) = 44843984081970031446.e-21;
      wi(11) = 48578046448352037528.e-21;
      wi(12) = 52070096091704461881.e-21;
      wi(13) = 55302735563728052549.e-21;
      wi(14) = 58259859877595495334.e-21;
      wi(15) = 60926736701561968039.e-21;
      wi(16) = 63290079733203854950.e-21;
      wi(17) = 65338114879181434984.e-21;
      wi(18) = 67060638906293652396.e-21;
      wi(19) = 68449070269366660985.e-21;
      wi(20) = 69496491861572578037.e-21;
      wi(21) = 70197685473558212587.e-21;
      wi(22) = 70549157789354068811.e-21;
      
    }

    else if (n==48) {
      xi(1) = -99877100725242611860.e-20;
      xi(2) = -99353017226635075755.e-20;
      xi(3) = -98412458372282685774.e-20;
      xi(4) = -97059159254624725046.e-20;
      xi(5) = -95298770316043086072.e-20;
      xi(6) = -93138669070655433311.e-20;
      xi(7) = -90587913671556967282.e-20;
      xi(8) = -87657202027424788591.e-20;
      xi(9) = -84358826162439353071.e-20;
      xi(10) = -80706620402944262708.e-20;
      xi(11) = -76715903251574033925.e-20;
      xi(12) = -72403413092381465467.e-20;
      xi(13) = -67787237963266390521.e-20;
      xi(14) = -62886739677651362400.e-20;
      xi(15) = -57722472608397270382.e-20;
      xi(16) = -52316097472223303368.e-20;
      xi(17) = -46690290475095840454.e-20;
      xi(18) = -40868648199071672992.e-20;
      xi(19) = -34875588629216073816.e-20;
      xi(20) = -28736248735545557674.e-20;
      xi(21) = -22476379039468906122.e-20;
      xi(22) = -16122235606889171806.e-20;
      xi(23) = -97004699209462698930.e-21;
      xi(24) = -32380170962869362033.e-21;
      wi(1) = 31533460523058386327.e-22;
      wi(2) = 73275539012762621024.e-22;
      wi(3) = 11477234579234539490.e-21;
      wi(4) = 15579315722943848728.e-21;
      wi(5) = 19616160457355527814.e-21;
      wi(6) = 23570760839324379141.e-21;
      wi(7) = 27426509708356948200.e-21;
      wi(8) = 31167227832798088902.e-21;
      wi(9) = 34777222564770438893.e-21;
      wi(10) = 38241351065830706317.e-21;
      wi(11) = 41545082943464749214.e-21;
      wi(12) = 44674560856694280419.e-21;
      wi(13) = 47616658492490474826.e-21;
      wi(14) = 50359035553854474958.e-21;
      wi(15) = 52890189485193667096.e-21;
      wi(16) = 55199503699984162868.e-21;
      wi(17) = 57277292100403215705.e-21;
      wi(18) = 59114839698395635746.e-21;
      wi(19) = 60704439165893880053.e-21;
      wi(20) = 62039423159892663904.e-21;
      wi(21) = 63114192286254025657.e-21;
      wi(22) = 63924238584648186624.e-21;
      wi(23) = 64466164435950082207.e-21;
      wi(24) = 64737696812683922503.e-21;
      
    }

    else if (n==52) {
      xi(1) = -99895111110395027809.e-20;
      xi(2) = -99447759092921602925.e-20;
      xi(3) = -98644619565154984065.e-20;
      xi(4) = -97488388422174450314.e-20;
      xi(5) = -95983182693308655253.e-20;
      xi(6) = -94134385364135905684.e-20;
      xi(7) = -91948612891642453989.e-20;
      xi(8) = -89433689053449532252.e-20;
      xi(9) = -86598616284606758524.e-20;
      xi(10) = -83453543232673453496.e-20;
      xi(11) = -80009728343046832433.e-20;
      xi(12) = -76279499519374496028.e-20;
      xi(13) = -72276209974998319368.e-20;
      xi(14) = -68014190422716770209.e-20;
      xi(15) = -63508697769524592430.e-20;
      xi(16) = -58775860497957906990.e-20;
      xi(17) = -53832620928582743838.e-20;
      xi(18) = -48696674569809607778.e-20;
      xi(19) = -43386406771876167031.e-20;
      xi(20) = -37920826911609366925.e-20;
      xi(21) = -32319500343480782550.e-20;
      xi(22) = -26602478360500182747.e-20;
      xi(23) = -20790226415636605969.e-20;
      xi(24) = -14903550860694918049.e-20;
      xi(25) = -89635244648900565489.e-21;
      xi(26) = -29914109797338766044.e-21;
      wi(1) = 26913169500471111189.e-22;
      wi(2) = 62555239629732768999.e-22;
      wi(3) = 98026345794627520617.e-22;
      wi(4) = 13315114982340960657.e-21;
      wi(5) = 16780023396300735678.e-21;
      wi(6) = 20184891507980792203.e-21;
      wi(7) = 23517513553984461590.e-21;
      wi(8) = 26765953746504013449.e-21;
      wi(9) = 29918581147143946641.e-21;
      wi(10) = 32964109089718797915.e-21;
      wi(11) = 35891634835097232942.e-21;
      wi(12) = 38690678310423978985.e-21;
      wi(13) = 41351219500560271679.e-21;
      wi(14) = 43863734259000407995.e-21;
      wi(15) = 46219228372784793508.e-21;
      wi(16) = 48409269744074896854.e-21;
      wi(17) = 50426018566342377218.e-21;
      wi(18) = 52262255383906993034.e-21;
      wi(19) = 53911406932757264751.e-21;
      wi(20) = 55367569669302652549.e-21;
      wi(21) = 56625530902368597191.e-21;
      wi(22) = 57680787452526827654.e-21;
      wi(23) = 58529561771813868550.e-21;
      wi(24) = 59168815466042970369.e-21;
      wi(25) = 59596260171248158258.e-21;
      wi(26) = 59810365745291860248.e-21;
      
    }

    else if (n==56) {
      xi(1) = -99909434380146558435.e-20;
      xi(2) = -99523122608106974722.e-20;
      xi(3) = -98829371554016151109.e-20;
      xi(4) = -97830170914025638338.e-20;
      xi(5) = -96528590190549018363.e-20;
      xi(6) = -94928647956196263565.e-20;
      xi(7) = -93035288024749630055.e-20;
      xi(8) = -90854362042065549085.e-20;
      xi(9) = -88392610832782754079.e-20;
      xi(10) = -85657643376274863540.e-20;
      xi(11) = -82657913214288165167.e-20;
      xi(12) = -79402692289386649803.e-20;
      xi(13) = -75902042270512890220.e-20;
      xi(14) = -72166783445018808352.e-20;
      xi(15) = -68208461269447045550.e-20;
      xi(16) = -64039310680700689427.e-20;
      xi(17) = -59672218277066332010.e-20;
      xi(18) = -55120682485553461875.e-20;
      xi(19) = -50398771838438171420.e-20;
      xi(20) = -45521081487845957895.e-20;
      xi(21) = -40502688092709127812.e-20;
      xi(22) = -35359103217495452097.e-20;
      xi(23) = -30106225386722066905.e-20;
      xi(24) = -24760290943433720397.e-20;
      xi(25) = -19337823863527525824.e-20;
      xi(26) = -13855584681037624201.e-20;
      xi(27) = -83305186822435374440.e-21;
      xi(28) = -27797035287275437094.e-21;
      wi(1) = 23238553757732154976.e-22;
      wi(2) = 54025222460153377612.e-22;
      wi(3) = 84690631633078876618.e-22;
      wi(4) = 11509824340383382175.e-21;
      wi(5) = 14515089278021471807.e-21;
      wi(6) = 17475512911400946507.e-21;
      wi(7) = 20381929882402572635.e-21;
      wi(8) = 23225351562565316937.e-21;
      wi(9) = 25996987058391952192.e-21;
      wi(10) = 28688268473822741730.e-21;
      wi(11) = 31290876747310447868.e-21;
      wi(12) = 33796767115611761295.e-21;
      wi(13) = 36198193872315186036.e-21;
      wi(14) = 38487734259247662487.e-21;
      wi(15) = 40658311384744517880.e-21;
      wi(16) = 42703216084667086511.e-21;
      wi(17) = 44616127652692283213.e-21;
      wi(18) = 46391133373001896762.e-21;
      wi(19) = 48022746793600258121.e-21;
      wi(20) = 49505924683047578920.e-21;
      wi(21) = 50836082617798480560.e-21;
      wi(22) = 52009109151741399843.e-21;
      wi(23) = 53021378524010763968.e-21;
      wi(24) = 53869761865714485709.e-21;
      wi(25) = 54551636870889421062.e-21;
      wi(26) = 55064895901762425796.e-21;
      wi(27) = 55407952503245123218.e-21;
      wi(28) = 55579746306514395846.e-21;
      
    }

    else if (n==60) {
      xi(1) = -99921012322743602204.e-20;
      xi(2) = -99584052511883817387.e-20;
      xi(3) = -98978789522222171737.e-20;
      xi(4) = -98106720175259818562.e-20;
      xi(5) = -96970178876505273373.e-20;
      xi(6) = -95572225583999610740.e-20;
      xi(7) = -93916627611642324950.e-20;
      xi(8) = -92007847617762755286.e-20;
      xi(9) = -89851031081004594194.e-20;
      xi(10) = -87451992264689831513.e-20;
      xi(11) = -84817198478592963249.e-20;
      xi(12) = -81953752616214575937.e-20;
      xi(13) = -78869373993226405457.e-20;
      xi(14) = -75572377530658568687.e-20;
      xi(15) = -72071651335573039944.e-20;
      xi(16) = -68376632738135543722.e-20;
      xi(17) = -64497282848947706781.e-20;
      xi(18) = -60444059704851036344.e-20;
      xi(19) = -56227890075394453918.e-20;
      xi(20) = -51860140005856974742.e-20;
      xi(21) = -47352584176170711111.e-20;
      xi(22) = -42717374158307838931.e-20;
      xi(23) = -37967005657679797715.e-20;
      xi(24) = -33114284826844819425.e-20;
      xi(25) = -28172293742326169169.e-20;
      xi(26) = -23154355137602933801.e-20;
      xi(27) = -18073996487342541724.e-20;
      xi(28) = -12944913539694500315.e-20;
      xi(29) = -77809333949536569419.e-21;
      xi(30) = -25959772301247798589.e-21;
      wi(1) = 20268119688737582862.e-22;
      wi(2) = 47127299269535689099.e-22;
      wi(3) = 73899311633454553761.e-22;
      wi(4) = 10047557182287984515.e-21;
      wi(5) = 12678166476815960077.e-21;
      wi(6) = 15274618596784799378.e-21;
      wi(7) = 17829901014207720191.e-21;
      wi(8) = 20337120729457286785.e-21;
      wi(9) = 22789516943997819862.e-21;
      wi(10) = 25180477621521248381.e-21;
      wi(11) = 27503556749924791635.e-21;
      wi(12) = 29752491500788945241.e-21;
      wi(13) = 31921219019296328949.e-21;
      wi(14) = 34003892724946422835.e-21;
      wi(15) = 35994898051084503067.e-21;
      wi(16) = 37888867569243444031.e-21;
      wi(17) = 39680695452380799470.e-21;
      wi(18) = 41365551235584755613.e-21;
      wi(19) = 42938892835935641954.e-21;
      wi(20) = 44396478795787113328.e-21;
      wi(21) = 45734379716114486647.e-21;
      wi(22) = 46948988848912204847.e-21;
      wi(23) = 48037031819971180964.e-21;
      wi(24) = 48995575455756835389.e-21;
      wi(25) = 49822035690550181011.e-21;
      wi(26) = 50514184532509374598.e-21;
      wi(27) = 51070156069855627405.e-21;
      wi(28) = 51488451500980933995.e-21;
      wi(29) = 51767943174910187544.e-21;
      wi(30) = 51907877631220639733.e-21;
      
    }

    else if (n==64) {
      xi(1) = -99930504173577213947.e-20;
      xi(2) = -99634011677195527930.e-20;
      xi(3) = -99101337147674432089.e-20;
      xi(4) = -98333625388462595699.e-20;
      xi(5) = -97332682778991096379.e-20;
      xi(6) = -96100879965205371896.e-20;
      xi(7) = -94641137485840281604.e-20;
      xi(8) = -92956917213193957580.e-20;
      xi(9) = -91052213707850280575.e-20;
      xi(10) = -88931544599511410585.e-20;
      xi(11) = -86599939815409281976.e-20;
      xi(12) = -84062929625258036275.e-20;
      xi(13) = -81326531512279755974.e-20;
      xi(14) = -78397235894334140761.e-20;
      xi(15) = -75281990726053189661.e-20;
      xi(16) = -71988185017161082685.e-20;
      xi(17) = -68523631305423324256.e-20;
      xi(18) = -64896547125465733986.e-20;
      xi(19) = -61115535517239325025.e-20;
      xi(20) = -57189564620263403428.e-20;
      xi(21) = -53127946401989454566.e-20;
      xi(22) = -48940314570705295748.e-20;
      xi(23) = -44636601725346408798.e-20;
      xi(24) = -40227015796399160370.e-20;
      xi(25) = -35722015833766811595.e-20;
      xi(26) = -31132287199021095616.e-20;
      xi(27) = -26468716220876741637.e-20;
      xi(28) = -21742364374000708415.e-20;
      xi(29) = -16964442042399281804.e-20;
      xi(30) = -12146281929612055447.e-20;
      xi(31) = -72993121787799039450.e-21;
      xi(32) = -24350292663424432509.e-21;
      wi(1) = 17832807216964305761.e-22;
      wi(2) = 41470332605624720784.e-22;
      wi(3) = 65044579689783580339.e-22;
      wi(4) = 88467598263639498816.e-22;
      wi(5) = 11168139460131126446.e-21;
      wi(6) = 13463047896718642240.e-21;
      wi(7) = 15726030476024719270.e-21;
      wi(8) = 17951715775697343343.e-21;
      wi(9) = 20134823153530209297.e-21;
      wi(10) = 22270173808383254264.e-21;
      wi(11) = 24352702568710873323.e-21;
      wi(12) = 26377469715054658680.e-21;
      wi(13) = 28339672614259483226.e-21;
      wi(14) = 30234657072402478867.e-21;
      wi(15) = 32057928354851553585.e-21;
      wi(16) = 33805161837141609392.e-21;
      wi(17) = 35472213256882383811.e-21;
      wi(18) = 37055128540240046040.e-21;
      wi(19) = 38550153178615629129.e-21;
      wi(20) = 39953741132720341387.e-21;
      wi(21) = 41262563242623528610.e-21;
      wi(22) = 42473515123653589007.e-21;
      wi(23) = 43583724529323453377.e-21;
      wi(24) = 44590558163756563060.e-21;
      wi(25) = 45491627927418144480.e-21;
      wi(26) = 46284796581314417296.e-21;
      wi(27) = 46968182816210017325.e-21;
      wi(28) = 47540165714830308662.e-21;
      wi(29) = 47999388596458307728.e-21;
      wi(30) = 48344762234802957170.e-21;
      wi(31) = 48575467441503426935.e-21;
      wi(32) = 48690957009139720383.e-21;
      
    }

    else if (n==72) {
      xi(1) = -99944993445296262425.e-20;
      xi(2) = -99710287164272906797.e-20;
      xi(3) = -99288495101680195780.e-20;
      xi(4) = -98680315237583044393.e-20;
      xi(5) = -97886877855723380415.e-20;
      xi(6) = -96909669799878047771.e-20;
      xi(7) = -95750524757769826743.e-20;
      xi(8) = -94411618527253794342.e-20;
      xi(9) = -92895464588091801987.e-20;
      xi(10) = -91204909268867146184.e-20;
      xi(11) = -89343126358809125431.e-20;
      xi(12) = -87313611129877890577.e-20;
      xi(13) = -85120173765443785094.e-20;
      xi(14) = -82766932202275431476.e-20;
      xi(15) = -80258304396929185201.e-20;
      xi(16) = -77599000029998256214.e-20;
      xi(17) = -74794011663283239067.e-20;
      xi(18) = -71848605366223496136.e-20;
      xi(19) = -68768310829046780749.e-20;
      xi(20) = -65558910981120105615.e-20;
      xi(21) = -62226431133946819101.e-20;
      xi(22) = -58777127669165559920.e-20;
      xi(23) = -55217476292771439740.e-20;
      xi(24) = -51554159877600352386.e-20;
      xi(25) = -47794055916894045718.e-20;
      xi(26) = -43944223612496073208.e-20;
      xi(27) = -40011890621916161270.e-20;
      xi(28) = -36004439489141923754.e-20;
      xi(29) = -31929393784671217133.e-20;
      xi(30) = -27794403980784760698.e-20;
      xi(31) = -23607233088575992506.e-20;
      xi(32) = -19375742083702605980.e-20;
      xi(33) = -15107875148221003033.e-20;
      xi(34) = -10811644756210281481.e-20;
      xi(35) = -64951166311857114075.e-21;
      xi(36) = -21663946035424044670.e-21;
      wi(1) = 14115163939753270721.e-22;
      wi(2) = 32831697746674940448.e-22;
      wi(3) = 51514360187894788637.e-22;
      wi(4) = 70102723218632486081.e-22;
      wi(5) = 88559960737053994464.e-22;
      wi(6) = 10685108165352501353.e-21;
      wi(7) = 12494165619873090060.e-21;
      wi(8) = 14279769054554032189.e-21;
      wi(9) = 16038564950285132556.e-21;
      wi(10) = 17767250789200705073.e-21;
      wi(11) = 19462580863294276956.e-21;
      wi(12) = 21121372216440556069.e-21;
      wi(13) = 22740510555035754038.e-21;
      wi(14) = 24316956064419165043.e-21;
      wi(15) = 25847749100655890089.e-21;
      wi(16) = 27330015738950934497.e-21;
      wi(17) = 28760973164701761061.e-21;
      wi(18) = 30137934895375479293.e-21;
      wi(19) = 31458315822561813978.e-21;
      wi(20) = 32719637064293846704.e-21;
      wi(21) = 33919530618286059497.e-21;
      wi(22) = 35055743807217870434.e-21;
      wi(23) = 36126143507637992986.e-21;
      wi(24) = 37128720154502899461.e-21;
      wi(25) = 38061591513802163834.e-21;
      wi(26) = 38923006216169663800.e-21;
      wi(27) = 39711347044834901782.e-21;
      wi(28) = 40425133971733970043.e-21;
      wi(29) = 41063026936075061102.e-21;
      wi(30) = 41623828360138598208.e-21;
      wi(31) = 42106485397586464147.e-21;
      wi(32) = 42510091910057720078.e-21;
      wi(33) = 42833890168338813667.e-21;
      wi(34) = 43077272274913699745.e-21;
      wi(35) = 43239781305222617485.e-21;
      wi(36) = 43321112165486537076.e-21;
      
    }

    else if (n==80) {
      xi(1) = -99955382265162826808.e-20;
      xi(2) = -99764986439820796834.e-20;
      xi(3) = -99422754096571174921.e-20;
      xi(4) = -98929130249972665638.e-20;
      xi(5) = -98284857273860426984.e-20;
      xi(6) = -97490914058572041817.e-20;
      xi(7) = -96548508904379327573.e-20;
      xi(8) = -95459076634363795066.e-20;
      xi(9) = -94224276130985166430.e-20;
      xi(10) = -92845987717243869507.e-20;
      xi(11) = -91326310257175780181.e-20;
      xi(12) = -89667557943877142823.e-20;
      xi(13) = -87872256767821427668.e-20;
      xi(14) = -85943140666311134922.e-20;
      xi(15) = -83883147358025523066.e-20;
      xi(16) = -81695413868146347481.e-20;
      xi(17) = -79383271750460544325.e-20;
      xi(18) = -76950242013504137461.e-20;
      xi(19) = -74400029758359727224.e-20;
      xi(20) = -71736518536209988023.e-20;
      xi(21) = -68963764434202760078.e-20;
      xi(22) = -66085989898611980174.e-20;
      xi(23) = -63107577304687196625.e-20;
      xi(24) = -60033062282975174315.e-20;
      xi(25) = -56867126812270978473.e-20;
      xi(26) = -53614592089713193202.e-20;
      xi(27) = -50280411188878498759.e-20;
      xi(28) = -46869661517054447704.e-20;
      xi(29) = -43387537083175609306.e-20;
      xi(30) = -39839340588196922702.e-20;
      xi(31) = -36230475349948731562.e-20;
      xi(32) = -32566437074770191462.e-20;
      xi(33) = -28852805488451185311.e-20;
      xi(34) = -25095235839227212049.e-20;
      xi(35) = -21299450285766613257.e-20;
      xi(36) = -17471229183264681256.e-20;
      xi(37) = -13616402280914388656.e-20;
      xi(38) = -97408398441584599063.e-21;
      xi(39) = -58504437152420668629.e-21;
      xi(40) = -19511383256793997654.e-21;
      wi(1) = 11449500037252453354.e-22;
      wi(2) = 26635335911449182005.e-22;
      wi(3) = 41803131248325406933.e-22;
      wi(4) = 56909224518929595161.e-22;
      wi(5) = 71929047688527908617.e-22;
      wi(6) = 86839452691609565673.e-22;
      wi(7) = 10161766041219293248.e-21;
      wi(8) = 11624114120416127748.e-21;
      wi(9) = 13068761592477329935.e-21;
      wi(10) = 14493508040524581339.e-21;
      wi(11) = 15896183583752925452.e-21;
      wi(12) = 17274652056269915690.e-21;
      wi(13) = 18626814208301626849.e-21;
      wi(14) = 19950610878140102324.e-21;
      wi(15) = 21244026115781505420.e-21;
      wi(16) = 22505090246332526056.e-21;
      wi(17) = 23731882865930076655.e-21;
      wi(18) = 24922535764115501728.e-21;
      wi(19) = 26075235767565117065.e-21;
      wi(20) = 27188227500486381444.e-21;
      wi(21) = 28259816057276862255.e-21;
      wi(22) = 29288369583267847694.e-21;
      wi(23) = 30272321759557980659.e-21;
      wi(24) = 31210174188114701643.e-21;
      wi(25) = 32100498673487773148.e-21;
      wi(26) = 32941939397645401383.e-21;
      wi(27) = 33733214984611522817.e-21;
      wi(28) = 34473120451753928794.e-21;
      wi(29) = 35160529044747593496.e-21;
      wi(30) = 35794393953416054603.e-21;
      wi(31) = 36373749905835978044.e-21;
      wi(32) = 36897714638276008839.e-21;
      wi(33) = 37365490238730490027.e-21;
      wi(34) = 37776364362001397490.e-21;
      wi(35) = 38129711314477638344.e-21;
      wi(36) = 38424993006959423185.e-21;
      wi(37) = 38661759774076463327.e-21;
      wi(38) = 38839651059051968932.e-21;
      wi(39) = 38958395962769531199.e-21;
      wi(40) = 39017813656306654811.e-21;
      
    }

    else if (n==88) {
      xi(1) = -99963083606662645367.e-20;
      xi(2) = -99805540804996391115.e-20;
      xi(3) = -99522317584752983826.e-20;
      xi(4) = -99113707813467295361.e-20;
      xi(5) = -98580218616336573201.e-20;
      xi(6) = -97922520312124670724.e-20;
      xi(7) = -97141441021462325102.e-20;
      xi(8) = -96237964604708261699.e-20;
      xi(9) = -95213229349011862112.e-20;
      xi(10) = -94068526328823504199.e-20;
      xi(11) = -92805297840849639853.e-20;
      xi(12) = -91425135510152027922.e-20;
      xi(13) = -89929778323076240540.e-20;
      xi(14) = -88321110405889090804.e-20;
      xi(15) = -86601158661395039865.e-20;
      xi(16) = -84772090209816777675.e-20;
      xi(17) = -82836209659047569960.e-20;
      xi(18) = -80795956199923475240.e-20;
      xi(19) = -78653900532828988584.e-20;
      xi(20) = -76412741628420688124.e-20;
      xi(21) = -74075303326860783021.e-20;
      xi(22) = -71644530779725048283.e-20;
      xi(23) = -69123486739088043741.e-20;
      xi(24) = -66515347698445317191.e-20;
      xi(25) = -63823399890332645027.e-20;
      xi(26) = -61051035145681898624.e-20;
      xi(27) = -58201746620128735758.e-20;
      xi(28) = -55279124392655547280.e-20;
      xi(29) = -52286850942114375373.e-20;
      xi(30) = -49228696507328672289.e-20;
      xi(31) = -46108514336619651061.e-20;
      xi(32) = -42930235832742437021.e-20;
      xi(33) = -39697865599349105025.e-20;
      xi(34) = -36415476395219828782.e-20;
      xi(35) = -33087204002619627511.e-20;
      xi(36) = -29717242016246430559.e-20;
      xi(37) = -26309836559336260023.e-20;
      xi(38) = -22869280933583131390.e-20;
      xi(39) = -19399910209614678858.e-20;
      xi(40) = -15906095764839421594.e-20;
      xi(41) = -12392239775547906165.e-20;
      xi(42) = -88627696702076058805.e-21;
      xi(43) = -53221325509403576522.e-21;
      xi(44) = -17747895902112098289.e-21;
      wi(1) = 94733513612529545899.e-23;
      wi(2) = 22040586988692204378.e-22;
      wi(3) = 34598684102906050894.e-22;
      wi(4) = 47114803663497454871.e-22;
      wi(5) = 59571849982300390278.e-22;
      wi(6) = 71953989055383030738.e-22;
      wi(7) = 84245466103046449651.e-22;
      wi(8) = 96430834163440844925.e-22;
      wi(9) = 10849469467224321456.e-21;
      wi(10) = 12042186598227368180.e-21;
      wi(11) = 13219730188112258690.e-21;
      wi(12) = 14380617607801932514.e-21;
      wi(13) = 15523385539474914063.e-21;
      wi(14) = 16646594210848201219.e-21;
      wi(15) = 17748828358570096473.e-21;
      wi(16) = 18828699173754781504.e-21;
      wi(17) = 19884846011190532102.e-21;
      wi(18) = 20915938130610818211.e-21;
      wi(19) = 21920676359974122618.e-21;
      wi(20) = 22897794734780677291.e-21;
      wi(21) = 23846062091859655017.e-21;
      wi(22) = 24764283620768778908.e-21;
      wi(23) = 25651302368961951718.e-21;
      wi(24) = 26506000699434738764.e-21;
      wi(25) = 27327301698855331284.e-21;
      wi(26) = 28114170534408613493.e-21;
      wi(27) = 28865615757635429249.e-21;
      wi(28) = 29580690553619349115.e-21;
      wi(29) = 30258493933943525335.e-21;
      wi(30) = 30898171871912197634.e-21;
      wi(31) = 31498918378604892320.e-21;
      wi(32) = 32059976518406388069.e-21;
      wi(33) = 32580639362732108686.e-21;
      wi(34) = 33060250880746700145.e-21;
      wi(35) = 33498206765953092528.e-21;
      wi(36) = 33893955197610259240.e-21;
      wi(37) = 34246997536020078737.e-21;
      wi(38) = 34556888950807084135.e-21;
      wi(39) = 34823238981399354993.e-21;
      wi(40) = 35045712029004261397.e-21;
      wi(41) = 35224027779459108533.e-21;
      wi(42) = 35357961556423843794.e-21;
      wi(43) = 35447344604470769706.e-21;
      wi(44) = 35492064301714545296.e-21;
      
    }

    else if (n==96) {
      xi(1) = -99968950458161677884.e-20;
      xi(2) = -99836436150026625442.e-20;
      xi(3) = -99598185401837168989.e-20;
      xi(4) = -99254388318724897532.e-20;
      xi(5) = -98805415112002250496.e-20;
      xi(6) = -98251725228772058939.e-20;
      xi(7) = -97593918962884178795.e-20;
      xi(8) = -96832679623798884729.e-20;
      xi(9) = -95968830381449102833.e-20;
      xi(10) = -95003272972082660508.e-20;
      xi(11) = -93937032620255700317.e-20;
      xi(12) = -92771245540840037901.e-20;
      xi(13) = -91507142091195428283.e-20;
      xi(14) = -90146063448296036815.e-20;
      xi(15) = -88689451757413437798.e-20;
      xi(16) = -87138850590440317660.e-20;
      xi(17) = -85495903346571022516.e-20;
      xi(18) = -83762351122975346770.e-20;
      xi(19) = -81940031074042640835.e-20;
      xi(20) = -80030874413847743608.e-20;
      xi(21) = -78036904386796125576.e-20;
      xi(22) = -75960234117661479314.e-20;
      xi(23) = -73803064374440921909.e-20;
      xi(24) = -71567681234896650925.e-20;
      xi(25) = -69256453664217173081.e-20;
      xi(26) = -66871831004391619210.e-20;
      xi(27) = -64416340378496711088.e-20;
      xi(28) = -61892584012546857021.e-20;
      xi(29) = -59303236477757208075.e-20;
      xi(30) = -56651041856139716841.e-20;
      xi(31) = -53938810832435743623.e-20;
      xi(32) = -51169417715466767359.e-20;
      xi(33) = -48345797392059635977.e-20;
      xi(34) = -45470942216774300864.e-20;
      xi(35) = -42547898840730054536.e-20;
      xi(36) = -39579764982890860329.e-20;
      xi(37) = -36569686147231363503.e-20;
      xi(38) = -33520852289262542262.e-20;
      xi(39) = -30436494435449635302.e-20;
      xi(40) = -27319881259104914149.e-20;
      xi(41) = -24174315616384001233.e-20;
      xi(42) = -21003131046056720360.e-20;
      xi(43) = -17809688236761860276.e-20;
      xi(44) = -14597371465489694199.e-20;
      xi(45) = -11369585011066592091.e-20;
      xi(46) = -81297495464425558994.e-21;
      xi(47) = -48812985136049731112.e-21;
      xi(48) = -16276744849602969579.e-21;
      wi(1) = 79647473413013308824.e-23;
      wi(2) = 18545712028943610772.e-22;
      wi(3) = 29096372771815468255.e-22;
      wi(4) = 39656312639089457047.e-22;
      wi(5) = 50134769158268190741.e-22;
      wi(6) = 60590278611871775766.e-22;
      wi(7) = 70960176225543354086.e-22;
      wi(8) = 81274503287260803750.e-22;
      wi(9) = 91484770246534778594.e-22;
      wi(10) = 10160502506419551876.e-21;
      wi(11) = 11162229814602627176.e-21;
      wi(12) = 12151660787661414925.e-21;
      wi(13) = 13128256444197094780.e-21;
      wi(14) = 14090951969013093016.e-21;
      wi(15) = 15038721680197214168.e-21;
      wi(16) = 15970562935363605962.e-21;
      wi(17) = 16885479435374326053.e-21;
      wi(18) = 17782502280253696814.e-21;
      wi(19) = 18660679624254468586.e-21;
      wi(20) = 19519081141356853600.e-21;
      wi(21) = 20356797152792231598.e-21;
      wi(22) = 21172939892395877365.e-21;
      wi(23) = 21966644438730298358.e-21;
      wi(24) = 22737069658335806966.e-21;
      wi(25) = 23483399085928249636.e-21;
      wi(26) = 24204841792364550405.e-21;
      wi(27) = 24900633222483566277.e-21;
      wi(28) = 25570036005349360697.e-21;
      wi(29) = 26212340735672413804.e-21;
      wi(30) = 26826866725591762076.e-21;
      wi(31) = 27412962726029242828.e-21;
      wi(32) = 27970007616848334438.e-21;
      wi(33) = 28497411065085385645.e-21;
      wi(34) = 28994614150555236543.e-21;
      wi(35) = 29461089958167905970.e-21;
      wi(36) = 29896344136328385984.e-21;
      wi(37) = 30299915420827593794.e-21;
      wi(38) = 30671376123669149014.e-21;
      wi(39) = 31010332586313837423.e-21;
      wi(40) = 31316425596861355813.e-21;
      wi(41) = 31589330770727168558.e-21;
      wi(42) = 31828758894411006535.e-21;
      wi(43) = 32034456231992663218.e-21;
      wi(44) = 32206204794030250669.e-21;
      wi(45) = 32343822568575928429.e-21;
      wi(46) = 32447163714064269364.e-21;
      wi(47) = 32516118713868835987.e-21;
      wi(48) = 32550614492363166242.e-21;
      
    }

    else if (n==104) {
      xi(1) = -99974025462451620190.e-20;
      xi(2) = -99861349070638500381.e-20;
      xi(3) = -99660327736123583184.e-20;
      xi(4) = -99361508866476523597.e-20;
      xi(5) = -98976332950492964529.e-20;
      xi(6) = -98504935190347884894.e-20;
      xi(7) = -97950469321598158611.e-20;
      xi(8) = -97301711615545710014.e-20;
      xi(9) = -96557677384453027357.e-20;
      xi(10) = -95730762749516840775.e-20;
      xi(11) = -94823341083979390414.e-20;
      xi(12) = -93824333288826003337.e-20;
      xi(13) = -92742344950081829375.e-20;
      xi(14) = -91576253747884907173.e-20;
      xi(15) = -90327594749081629450.e-20;
      xi(16) = -88997023646288045179.e-20;
      xi(17) = -87586186756351095122.e-20;
      xi(18) = -86096161127797694571.e-20;
      xi(19) = -84528340899234132869.e-20;
      xi(20) = -82884124840218026839.e-20;
      xi(21) = -81165006660051550143.e-20;
      xi(22) = -79372537644043611522.e-20;
      xi(23) = -77508338027046977385.e-20;
      xi(24) = -75574092471832298179.e-20;
      xi(25) = -73571549017801969428.e-20;
      xi(26) = -71502517397363087234.e-20;
      xi(27) = -69368867439397065602.e-20;
      xi(28) = -67172527368333210088.e-20;
      xi(29) = -64915482063111896041.e-20;
      xi(30) = -62599771263251515327.e-20;
      xi(31) = -60227487725540047159.e-20;
      xi(32) = -57800775332757748445.e-20;
      xi(33) = -55321827156203442307.e-20;
      xi(34) = -52792883473767725346.e-20;
      xi(35) = -50216229745345025782.e-20;
      xi(36) = -47594194547413982495.e-20;
      xi(37) = -44929147468652664376.e-20;
      xi(38) = -42223496968490362523.e-20;
      xi(39) = -39479688200531193751.e-20;
      xi(40) = -36700200802816507262.e-20;
      xi(41) = -33887546656923060512.e-20;
      xi(42) = -31044267617922098404.e-20;
      xi(43) = -28172933217250806882.e-20;
      xi(44) = -25276138340572094235.e-20;
      xi(45) = -22356500882721258998.e-20;
      xi(46) = -19416659381858811954.e-20;
      xi(47) = -16459270634967512819.e-20;
      xi(48) = -13487007296848542737.e-20;
      xi(49) = -10502555464786646703.e-20;
      xi(50) = -75086122510670317989.e-21;
      xi(51) = -45078833455377862647.e-21;
      xi(52) = -15030805704205808070.e-21;
      wi(1) = 80537301359223283137.e-23;
      wi(2) = 88612279172572261915.e-23;
      wi(3) = 29505344364932164953.e-22;
      wi(4) = 29860447460835539614.e-22;
      wi(5) = 44717827658964106576.e-22;
      wi(6) = 54373252304872340492.e-22;
      wi(7) = 63840289565270698417.e-22;
      wi(8) = 64909501207992000200.e-22;
      wi(9) = 83745523985271410846.e-22;
      wi(10) = 85549571290963629793.e-22;
      wi(11) = 94157585502386417660.e-22;
      wi(12) = 10364141417799710457.e-21;
      wi(13) = 11288156174522055729.e-21;
      wi(14) = 12085767303575836382.e-21;
      wi(15) = 12898236413125625150.e-21;
      wi(16) = 13709097817130044944.e-21;
      wi(17) = 14507175337544532779.e-21;
      wi(18) = 15291431851904123958.e-21;
      wi(19) = 16062483287466054902.e-21;
      wi(20) = 16819201771329851973.e-21;
      wi(21) = 17560584847393953123.e-21;
      wi(22) = 18286096131447045103.e-21;
      wi(23) = 18995087062019415460.e-21;
      wi(24) = 19686910416592402212.e-21;
      wi(25) = 20360942228950003547.e-21;
      wi(26) = 21016573507412752623.e-21;
      wi(27) = 21653211654050478214.e-21;
      wi(28) = 22270281336456832151.e-21;
      wi(29) = 22867224894903991760.e-21;
      wi(30) = 23443502859080308868.e-21;
      wi(31) = 23998594434228088202.e-21;
      wi(32) = 24531997972222651431.e-21;
      wi(33) = 25043231424904261260.e-21;
      wi(34) = 25531832779704979249.e-21;
      wi(35) = 25997360477175197755.e-21;
      wi(36) = 26439393810028101126.e-21;
      wi(37) = 26857533303339869813.e-21;
      wi(38) = 27251401075562153020.e-21;
      wi(39) = 27620641180020444979.e-21;
      wi(40) = 27964919926589683664.e-21;
      wi(41) = 28283926183256317188.e-21;
      wi(42) = 28577371657294272064.e-21;
      wi(43) = 28844991155800689526.e-21;
      wi(44) = 29086542825355955653.e-21;
      wi(45) = 29301808370591421861.e-21;
      wi(46) = 29490593251467277733.e-21;
      wi(47) = 29652726859082281213.e-21;
      wi(48) = 29788062669856454763.e-21;
      wi(49) = 29896478377947402691.e-21;
      wi(50) = 29977876005780577121.e-21;
      wi(51) = 30032181992593600121.e-21;
      wi(52) = 30059347260914619701.e-21;
      
    }

    else if (n==112) {
      xi(1) = -80488595332676199779.e-20;
      xi(2) = -10160478627668805531.e-19;
      xi(3) = -95439318547860299986.e-20;
      xi(4) = -91546298836742760945.e-20;
      xi(5) = -10340513562261696462.e-19;
      xi(6) = -92662585778661689507.e-20;
      xi(7) = -90489298887912825445.e-20;
      xi(8) = -10131949463069062515.e-19;
      xi(9) = -90464129755666057124.e-20;
      xi(10) = -86603614404095298415.e-20;
      xi(11) = -91636186250656082289.e-20;
      xi(12) = -10121834958204258714.e-19;
      xi(13) = -10078277129241111043.e-19;
      xi(14) = -87962237189807963118.e-20;
      xi(15) = -91677434744694723682.e-20;
      xi(16) = -90487658323691933037.e-20;
      xi(17) = -89271920604983483215.e-20;
      xi(18) = -87962151397116941034.e-20;
      xi(19) = -86599374258026102045.e-20;
      xi(20) = -85171753748871651951.e-20;
      xi(21) = -83675665375815917519.e-20;
      xi(22) = -82114054838778695246.e-20;
      xi(23) = -80488575638784496812.e-20;
      xi(24) = -78800297986553099992.e-20;
      xi(25) = -77050563606913539332.e-20;
      xi(26) = -75240747787024521741.e-20;
      xi(27) = -73372261848570062821.e-20;
      xi(28) = -71446562385808640412.e-20;
      xi(29) = -69465151030101343512.e-20;
      xi(30) = -67429572842563396177.e-20;
      xi(31) = -65341415102530449682.e-20;
      xi(32) = -63202306093043949251.e-20;
      xi(33) = -61013913826811007054.e-20;
      xi(34) = -58777944746249047206.e-20;
      xi(35) = -56496142392748458796.e-20;
      xi(36) = -54170286047113786527.e-20;
      xi(37) = -51802189342131382164.e-20;
      xi(38) = -49393698848352747124.e-20;
      xi(39) = -46946692634193765874.e-20;
      xi(40) = -44463078801473008955.e-20;
      xi(41) = -41944793997530961580.e-20;
      xi(42) = -39393801905090372790.e-20;
      xi(43) = -36812091711035275365.e-20;
      xi(44) = -34201676555302667702.e-20;
      xi(45) = -31564591961096358173.e-20;
      xi(46) = -28902894247647038221.e-20;
      xi(47) = -26218658926756261480.e-20;
      xi(48) = -23513979084374651836.e-20;
      xi(49) = -20790963748476333913.e-20;
      xi(50) = -18051736244502265770.e-20;
      xi(51) = -15298432539654847416.e-20;
      xi(52) = -12533199577334872554.e-20;
      xi(53) = -97581936030195779185.e-21;
      xi(54) = -69755784828872187578.e-21;
      xi(55) = -41875240164992552336.e-21;
      xi(56) = -13962042448558683275.e-21;
      wi(1) = 16569636576973647761.e-21;
      wi(2) = 14402201696983724428.e-39;
      wi(3) = 34709520255747461200.e-22;
      wi(4) = 56778376879010079815.e-21;
      wi(5) = 37731161063754411165.e-47;
      wi(6) = 45632271760860013289.e-22;
      wi(7) = 11520495792117312626.e-21;
      wi(8) = 49246746812766292842.e-38;
      wi(9) = 14372045777313838419.e-21;
      wi(10) = 13870830688358627358.e-21;
      wi(11) = -37371950806309312086.e-21;
      wi(12) = 18775853756740260905.e-37;
      wi(13) = 12261054337390208321.e-34;
      wi(14) = 13637759622784054281.e-21;
      wi(15) = 20686815045027872828.e-21;
      wi(16) = 84327272746907118565.e-22;
      wi(17) = 11759939136603779624.e-21;
      wi(18) = 13945650387627380053.e-21;
      wi(19) = 14124057734648147892.e-21;
      wi(20) = 14634623094715024217.e-21;
      wi(21) = 15300784134556015037.e-21;
      wi(22) = 15940621199485025417.e-21;
      wi(23) = 16571845899450453236.e-21;
      wi(24) = 17192076898900235720.e-21;
      wi(25) = 17800021943506593268.e-21;
      wi(26) = 18393902644953860296.e-21;
      wi(27) = 18973392242557999400.e-21;
      wi(28) = 19538093130310137860.e-21;
      wi(29) = 20087558367323984780.e-21;
      wi(30) = 20621359698898228923.e-21;
      wi(31) = 21139081089410710408.e-21;
      wi(32) = 21640318863160403311.e-21;
      wi(33) = 22124682168949492325.e-21;
      wi(34) = 22591793313385707526.e-21;
      wi(35) = 23041288057354547180.e-21;
      wi(36) = 23472815898299421860.e-21;
      wi(37) = 23886040343752650004.e-21;
      wi(38) = 24280639173689438499.e-21;
      wi(39) = 24656304691784812879.e-21;
      wi(40) = 25012743965344867773.e-21;
      wi(41) = 25349679053726235917.e-21;
      wi(42) = 25666847225065297520.e-21;
      wi(43) = 25964001161148146161.e-21;
      wi(44) = 26240909150261532990.e-21;
      wi(45) = 26497355267874394909.e-21;
      wi(46) = 26733139545009064708.e-21;
      wi(47) = 26948078124170862883.e-21;
      wi(48) = 27142003402714474402.e-21;
      wi(49) = 27314764163535311670.e-21;
      wi(50) = 27466225692983949724.e-21;
      wi(51) = 27596269885911683746.e-21;
      wi(52) = 27704795337765294488.e-21;
      wi(53) = 27791717423659206498.e-21;
      wi(54) = 27856968364363379146.e-21;
      wi(55) = 27900497279155473642.e-21;
      wi(56) = 27922270225496082396.e-21;

        } else {
            warn("case n= %d is not provided in GauLeg", n);
            return number + 1;
        }
#endif // 1
        for (int i = n; i >= (n + 3)/2; --i) {
            xi(i) = -xi(n + 1 - i);
            wi(i) =  wi(n + 1 - i);
        } // i
        return n;
    } // Gauss_Legendre_quadrature


    int Gauss_Fermi_Dirac_quadrature(double x[], double w[], unsigned const n, int const echo) {
#if 1
    if (n==0) {
        return 0;
    } else
    if (n==1) {
      xi(1) = -49817229548128141768.e-20;
      wi(1) = 10000000000000031192.e-19;
      

        } else

    if (n==2) {
      xi(1) = -78465071850839016234.e-20;
      xi(2) = -20091536266094051757.e-20;
      wi(1) = 50923235990870048433.e-20;
      wi(2) = 49076764009130263488.e-20;
      

        } else

    if (n==3) {
      xi(1) = -88288518955458358024.e-20;
      xi(2) = -48117621892777473749.e-20;
      xi(3) = -88198184413497647625.e-21;
      wi(1) = 28858444436509900908.e-20;
      wi(2) = 45966895698954759346.e-20;
      wi(3) = 25174659864535651667.e-20;
      

        } else

    if (n==4) {
      xi(1) = -92613063531202843773.e-20;
      xi(2) = -64918327008663578157.e-20;
      xi(3) = -28982568853420020298.e-20;
      xi(4) = -24595209663255169680.e-21;
      wi(1) = 18501429405165520392.e-20;
      wi(2) = 34614391006511784214.e-20;
      wi(3) = 34152482191988153127.e-20;
      wi(4) = 12731697396334854188.e-20;
      

        } else

    if (n==5) {
      xi(1) = -94875333872503463082.e-20;
      xi(2) = -74805843506753178608.e-20;
      xi(3) = -45504655263391074765.e-20;
      xi(4) = -16657582360358973599.e-20;
      xi(5) = 27402283545708211900.e-21;
      wi(1) = 12939804504572789754.e-20;
      wi(2) = 26102400189213290231.e-20;
      wi(3) = 30851911091450589451.e-20;
      wi(4) = 24746815229701880449.e-20;
      wi(5) = 53590689850617620359.e-21;
      

        } else

    if (n==6) {
      xi(1) = -96204950250095729781.e-20;
      xi(2) = -80971428101130972258.e-20;
      xi(3) = -57293627456482418171.e-20;
      xi(4) = -30755197635518367504.e-20;
      xi(5) = -82123839469384988331.e-21;
      xi(6) = 83748358371240941581.e-21;
      wi(1) = 96268650841705383829.e-21;
      wi(2) = 20246201047059595265.e-20;
      wi(3) = 26160719441051813381.e-20;
      wi(4) = 25781980698475975536.e-20;
      wi(5) = 16683001513553609336.e-20;
      wi(6) = 15012322156887800205.e-21;
      

        } else

    if (n==7) {
      xi(1) = -97053934379083423143.e-20;
      xi(2) = -85045695849615413757.e-20;
      xi(3) = -65665104053460540522.e-20;
      xi(4) = -42357896269371657364.e-20;
      xi(5) = -19472732441816555564.e-20;
      xi(6) = -19669621223691542539.e-21;
      xi(7) = 15142830586888806919.e-20;
      wi(1) = 74948008822570509041.e-21;
      wi(2) = 16170863905729061704.e-20;
      wi(3) = 22007120289205973485.e-20;
      wi(4) = 23880411919774885893.e-20;
      wi(5) = 20952460047488907594.e-20;
      wi(6) = 92465405554445737538.e-21;
      wi(7) = 24780240009985858690.e-22;
      

        } else

    if (n==8) {
      xi(1) = -97630544447925725992.e-20;
      xi(2) = -87873822716479965943.e-20;
      xi(3) = -71736329217593360204.e-20;
      xi(4) = -51463306578144813387.e-20;
      xi(5) = -29967081434747298359.e-20;
      xi(6) = -10763455942936048359.e-20;
      xi(7) = 35963113675701677498.e-21;
      xi(8) = 23003149140664609750.e-20;
      wi(1) = 60394634019629989770.e-21;
      wi(2) = 13252509350880929004.e-20;
      wi(3) = 18643612522057003210.e-20;
      wi(4) = 21413715867867937533.e-20;
      wi(5) = 21005092708864293339.e-20;
      wi(6) = 16003068683842947897.e-20;
      wi(7) = 36159126989806650464.e-21;
      wi(8) = 26624765543536915040.e-23;
      

        } else

    if (n==9) {
      xi(1) = -98041275487012188695.e-20;
      xi(2) = -89918326179154863440.e-20;
      xi(3) = -76254129548477842110.e-20;
      xi(4) = -58579104527384144901.e-20;
      xi(5) = -38924212142470946276.e-20;
      xi(6) = -19724340764961096691.e-20;
      xi(7) = -40039281758884590381.e-21;
      xi(8) = 97228170103579374416.e-21;
      xi(9) = 31678885353558278864.e-20;
      wi(1) = 49992516372028853833.e-21;
      wi(2) = 11099301824870447793.e-20;
      wi(3) = 15971411690431220541.e-20;
      wi(4) = 19037877203046567198.e-20;
      wi(5) = 19869087157813151863.e-20;
      wi(6) = 17972334325952047726.e-20;
      wi(7) = 10203571121909080322.e-20;
      wi(8) = 84501828581921130722.e-22;
      wi(9) = 21467529556997868476.e-24;
      

        } else

    if (n==10) {
      xi(1) = -98345122025502045873.e-20;
      xi(2) = -91446749996879318119.e-20;
      xi(3) = -79700500547314513626.e-20;
      xi(4) = -64189534981349313375.e-20;
      xi(5) = -46376588343242516012.e-20;
      xi(6) = -28030431525349494354.e-20;
      xi(7) = -11327091328726333942.e-20;
      xi(8) = 17437648086722052805.e-21;
      xi(9) = 16877498338102917782.e-20;
      xi(10) = 40960465258252015313.e-20;
      wi(1) = 42278597323639457484.e-21;
      wi(2) = 94666349251635366832.e-21;
      wi(3) = 13843777024241956101.e-20;
      wi(4) = 16932936699837666261.e-20;
      wi(5) = 18398357022114735352.e-20;
      wi(6) = 17939886390638648260.e-20;
      wi(7) = 14468854182396060463.e-20;
      wi(8) = 46026485095922891703.e-21;
      wi(9) = 11890402956686871419.e-22;
      wi(10) = 14148408460516817666.e-25;
      

        } else

    if (n==11) {
      xi(1) = -98576901837451635280.e-20;
      xi(2) = -92621727156102677473.e-20;
      xi(3) = -82389243156123939088.e-20;
      xi(4) = -68670708816882492198.e-20;
      xi(5) = -52549052940365991088.e-20;
      xi(6) = -35349156561982307316.e-20;
      xi(7) = -18652071146560858606.e-20;
      xi(8) = -45389164233559550280.e-21;
      xi(9) = 76984180593432347734.e-21;
      xi(10) = 24899533750455431614.e-20;
      xi(11) = 50711636785486806957.e-20;
      wi(1) = 36383684790132198923.e-21;
      wi(2) = 81985364434128201418.e-21;
      wi(3) = 12133566247788805356.e-20;
      wi(4) = 15122112006362489825.e-20;
      wi(5) = 16900090791849557413.e-20;
      wi(6) = 17240157268363508589.e-20;
      wi(7) = 15745585899461757802.e-20;
      wi(8) = 97600157144810676257.e-21;
      wi(9) = 12496828256639735424.e-21;
      wi(10) = 11876318920871395759.e-23;
      wi(11) = 80046822403386311030.e-27;
      

        } else

    if (n==12) {
      xi(1) = -98758247347129831371.e-20;
      xi(2) = -93546465146779806654.e-20;
      xi(3) = -84528996754470930223.e-20;
      xi(4) = -72299594230844519839.e-20;
      xi(5) = -57679398168141327066.e-20;
      xi(6) = -41683730779892996801.e-20;
      xi(7) = -25514627335790291149.e-20;
      xi(8) = -10710838211747769681.e-20;
      xi(9) = 12720145729326415607.e-21;
      xi(10) = 14540842218988328389.e-20;
      xi(11) = 33552500235752414908.e-20;
      xi(12) = 60838109964484063119.e-20;
      wi(1) = 31765161579790701148.e-21;
      wi(2) = 71927618746964313778.e-21;
      wi(3) = 10742555378156694842.e-20;
      wi(4) = 13578811351554214795.e-20;
      wi(5) = 15492042553417744038.e-20;
      wi(6) = 16300300254834219520.e-20;
      wi(7) = 15784577013790806216.e-20;
      wi(8) = 12921482926208917372.e-20;
      wi(9) = 46096943233133302568.e-21;
      wi(10) = 20030610755774790850.e-22;
      wi(11) = 95165705752725893549.e-25;
      wi(12) = 40143360822128708729.e-28;
      

        } else

    if (n==13) {
      xi(1) = -98903182721370020265.e-20;
      xi(2) = -94288936524363459773.e-20;
      xi(3) = -86261843870640242196.e-20;
      xi(4) = -75277808759167753869.e-20;
      xi(5) = -61972590294795871779.e-20;
      xi(6) = -47139332563986024748.e-20;
      xi(7) = -31718188942187627557.e-20;
      xi(8) = -16854863011308355787.e-20;
      xi(9) = -41195843159851553906.e-21;
      xi(10) = 71957380142115164738.e-21;
      xi(11) = 22223926926874000328.e-20;
      xi(12) = 42682885634093164862.e-20;
      xi(13) = 71270930856714354732.e-20;
      wi(1) = 28069991026027589482.e-21;
      wi(2) = 63803895087070663653.e-21;
      wi(3) = 95973484361405430270.e-21;
      wi(4) = 12264378189747678145.e-20;
      wi(5) = 14213612346123977130.e-20;
      wi(6) = 15296686007570952707.e-20;
      wi(7) = 15358437552921000921.e-20;
      wi(8) = 14007635729175637795.e-20;
      wi(9) = 87531230524252970103.e-21;
      wi(10) = 12989730151883234012.e-21;
      wi(11) = 22351943999969127535.e-23;
      wi(12) = 65097139765619073344.e-26;
      wi(13) = 18257341724040876662.e-29;
      

        } else

    if (n==14) {
      xi(1) = -99021130855943209687.e-20;
      xi(2) = -94895368426058288869.e-20;
      xi(3) = -87686856465753704289.e-20;
      xi(4) = -77752669471002194917.e-20;
      xi(5) = -65594116901081876554.e-20;
      xi(6) = -51841232227159879604.e-20;
      xi(7) = -37243750660439082187.e-20;
      xi(8) = -22693429290756856295.e-20;
      xi(9) = -93940943648510570987.e-21;
      xi(10) = 16521198218716065629.e-21;
      xi(11) = 13919799114797561344.e-20;
      xi(12) = 30521886852802066309.e-20;
      xi(13) = 52192337126752562221.e-20;
      xi(14) = 81957965081548293179.e-20;
      wi(1) = 25060310888021301605.e-21;
      wi(2) = 57137272611562033779.e-21;
      wi(3) = 86434450014324433897.e-21;
      wi(4) = 11141118228632175288.e-20;
      wi(5) = 13070790263291078499.e-20;
      wi(6) = 14310195071194851995.e-20;
      wi(7) = 14737968606274298328.e-20;
      wi(8) = 14154903694980505066.e-20;
      wi(9) = 11456160782223814050.e-20;
      wi(10) = 40466499493397342820.e-21;
      wi(11) = 21701008894932486895.e-22;
      wi(12) = 19960253076851250807.e-24;
      wi(13) = 39376501060604877095.e-27;
      wi(14) = 76596142918862399780.e-31;
      

        } else

    if (n==15) {
      xi(1) = -99118619138431485634.e-20;
      xi(2) = -95398089203095832045.e-20;
      xi(3) = -88874665207045485764.e-20;
      xi(4) = -79832886799647722652.e-20;
      xi(5) = -68674462209286747178.e-20;
      xi(6) = -55907326778454372362.e-20;
      xi(7) = -42138595122137487519.e-20;
      xi(8) = -28083407355763995168.e-20;
      xi(9) = -14649293944496725019.e-20;
      xi(10) = -30865949117072113052.e-21;
      xi(11) = 75989566859912966734.e-21;
      xi(12) = 21425891814116860148.e-20;
      xi(13) = 39280262275215780450.e-20;
      xi(14) = 62012182191671475949.e-20;
      xi(15) = 92858877219218103945.e-20;
      wi(1) = 22570991165870390473.e-21;
      wi(2) = 51589746641923392000.e-21;
      wi(3) = 78401918844466166239.e-21;
      wi(4) = 10176234626640128024.e-20;
      wi(5) = 12055819130110177262.e-20;
      wi(6) = 13377324647273569326.e-20;
      wi(7) = 14041818603247829422.e-20;
      wi(8) = 13919569003129657925.e-20;
      wi(9) = 12562361445602688222.e-20;
      wi(10) = 74852662340708470150.e-21;
      wi(11) = 10996744175647251144.e-21;
      wi(12) = 25513307315040157893.e-23;
      wi(13) = 15270418102934789627.e-25;
      wi(14) = 21560859319293022163.e-28;
      wi(15) = 30032040385910287756.e-32;
      

        } else

    if (n==16) {
      xi(1) = -99200289748411473927.e-20;
      xi(2) = -95820266378296634182.e-20;
      xi(3) = -89876661129475763142.e-20;
      xi(4) = -81599671254865832971.e-20;
      xi(5) = -71315812647978444249.e-20;
      xi(6) = -59440032425488487666.e-20;
      xi(7) = -46470396871945791541.e-20;
      xi(8) = -32991653294098863600.e-20;
      xi(9) = -19716091326862980561.e-20;
      xi(10) = -76605243443508959615.e-21;
      xi(11) = 26155046503992069925.e-21;
      xi(12) = 14307776307824938682.e-20;
      xi(13) = 29506185654032182160.e-20;
      xi(14) = 48403577800553841578.e-20;
      xi(15) = 72091584865612160132.e-20;
      xi(16) = 10394188783181811718.e-19;
      wi(1) = 20484388078614008045.e-21;
      wi(2) = 46916532350372347409.e-21;
      wi(3) = 71569877291069983495.e-21;
      wi(4) = 93424466379672137196.e-21;
      wi(5) = 11156011364306951297.e-20;
      wi(6) = 12512553084306063601.e-20;
      wi(7) = 13329704953113185969.e-20;
      wi(8) = 13510959073859290681.e-20;
      wi(9) = 12840858805365359846.e-20;
      wi(10) = 10016528657871746742.e-20;
      wi(11) = 32102655847303900301.e-21;
      wi(12) = 18115418480524121495.e-22;
      wi(13) = 24274994772381143993.e-24;
      wi(14) = 10371321943363515335.e-26;
      wi(15) = 10868941709467004901.e-29;
      wi(16) = 11117372791599461059.e-33;
        } else {
            warn("number must be in [0, 16], found= %d", n);
            return n + 1;
        }
#endif // 1

#undef xi
#undef wi
       return n;
    } // Gauss_Fermi_Dirac_quadrature 



#ifdef    NO_UNIT_TESTS
    status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

    status_t test_integrator(int const echo=3) {
        if (echo > 1) std::printf("\n#\n# %s\n", __func__);
        status_t stat(0);
        { // scope
            // spectrum of an isolated box  0.648721 1.29538_x3 1.94203_x3 2.34449_x3 2.58869 2.99114_x6 3.55103_x3 3.6378_x3 ...
            //  ... 4.04025_x3 4.19768_x6 4.68691_x3 4.84434_x3 5.24679_x6 5.73601 5.89345_x6 6.45333_x3 6.94256_x3 7.09999_x3 8.1491_x3 9.35564

            double E_Fermi{1.0};
            std::vector<double> xyzZinso(0); // no atoms
            real_space::grid_t gc(4, 4, 4); // one block, isolated BCs by default, grid spacing 1.0
            parallel_poisson::load_balancing_t const lb(gc, MPI_COMM_WORLD, 4, echo);
            parallel_poisson::parallel_grid_t const pg(gc, lb, echo, "Interpolation");
            view2D<double> V_coarse(pg.n_local(), 4*4*4, 0.5);
            view2D<double> rhov_new(pg.n_local(), 8*8*8, 0.0);
            std::vector<uint32_t> num(0);
            data_list<double> atom_mat(num);
            std::vector<int32_t> numax_prj(0, 0);
            std::vector<double> sigma_prj(0, 1.);
            Integrator integrator(gc, xyzZinso, echo);
            if (echo > 1) std::printf("# %s: Integrator constructed\n\n", __func__);
            stat += integrator.integrate(rhov_new[0], E_Fermi, V_coarse[0], atom_mat, numax_prj, sigma_prj, pg, 1., 1., echo);
        } // scope (so all destructors belonging to this test are called before the next log message)
        if (echo > 1) std::printf("# %s: Integrator.integrate executed\n\n", __func__);
        return stat;
    } // test_integrator

    status_t test_energy_mesh(int const echo=5) {
        double c_ref{0};
        { // scope: compute reference value by integrating over real axis [-1, 0]
            int constexpr m = 1000;
            auto const dx = 1.0/m;
            for (int i{0}; i < 1000; ++i) {
                c_ref += std::cos(-(i + .5)*dx)*dx;
            } // i
        } // scope
        auto const reference = c_ref;

        auto const nBot = std::min(unsigned(control::get("energy_contour.bottom",    9.)), 112u);
        auto const nPar = std::min(unsigned(control::get("energy_contour.parallel", 33.)), 112u);
        auto const nFer = std::min(unsigned(control::get("energy_contour.fermidirac", 9.)), 16u);
        auto const nPol =               int(control::get("energy_contour.matsubara", 1.));
        double const eBot = -1, kBT = 1.008034e-2; // == 2e4 Kelvin
        for (unsigned nb{0}; nb <= nBot; ++nb) {
        for (unsigned np{0}; np <= nPar; ++np) {
        for (unsigned nf{0}; nf <= nFer; ++nf) {
        for (int nm{1}; nm <= nPol; ++nm) {
            std::vector<Complex> energy_weights;
            auto const energy_mesh = get_energy_mesh(energy_weights, kBT, eBot, nb, np, nf, nm, echo/4);
            if (echo > 1) std::printf("# energy mesh with %ld points generated\n", energy_mesh.size());
            // integrate a simple but holomorphic function of which the integral over the real axis [-1, 0] is known
            auto const nE = energy_mesh.size();
            assert(nE == energy_weights.size());
            Complex c(0);
            for (size_t iE{0}; iE < nE; ++iE) {
                auto const x = energy_mesh.at(iE), wgt = energy_weights.at(iE);
                c += std::cos(x)*wgt;
            } // iE
            if (echo > 0) std::printf("# %s(nBot=%d, nPar=%d, nFer=%d, nPol=%d) Integral(cos)= (%g, %g), reference= %g\n",
                                         __func__, nb, np, nf, nm, c.real(), c.imag(), reference);
        }}}} // nb np nf nm
        return 0;
    } // test_energy_mesh

    view2D<double> Legendre_polynomials(unsigned const n, int const echo=0) {
        view2D<double> a(n, n, 0.0);
        if (n > 0) { set(a[0], n, 0.0); a(0,0) = 1; }
        if (n > 1) { set(a[1], n, 0.0); a(1,1) = 1; }
        // (n+1)*P_{n+1}(x) = (2*n+1)*x*P_{n}(x) - n*P_{n-1}(x)
        for (int k{1}; k < n - 1; ++k) {
            set(a[k + 1], n, a[k - 1], -double(k));
            add_product(&a(k + 1,1), n - 1, a[k], 2*k + 1.);
            scale(a[k + 1], n, 1./(n + 1));
        } // k
        if (echo > 9) {
            for (int k{0}; k < std::min(n, 10u); ++k) {
                std::printf("# %s n= %d: \t", __func__, k); printf_vector(" %g", a[k], k + 1);
            } // k
        } // echo
        return a;
    } // Legendre_polynomials

    status_t test_Gauss_Legendre_quadrature(int const echo=5) {
        status_t stat(0);
        double x[120], w[120];
        auto const polys = Legendre_polynomials(120, echo);
        for (int n{0}; n <= 12; ++n) {
            auto const nn = Gauss_Legendre_quadrature(x, w, n, echo*(n < 33));
            if (nn == n) {
                // integrate the Legendre polynomials using the quadrature weights
                std::vector<double> s(n, 0.0);
                for (int i{0}; i < n; ++i) { // for each point in the quadrature
                    double xipow{1};
                    for (int k{0}; k < n; ++k) { // for each polynomial degree
                        for (int j{0}; j <= k; ++j) { // for each polynomial
                            s[j] += polys(k,j)*xipow*w[i];
                        } // j
                        xipow *= x[i]; // prepare (x[i])^k for k := (k+1)
                    } // k
                } // i
                std::printf("# integrals for n=%d ", n); printf_vector(" %.6f", s);
            } // success
        } // n
        return stat;
    } // test_Gauss_Legendre_quadrature

    status_t all_tests(int const echo) {
        int const which = control::get("energy_contour.select.test", -1.);
        status_t stat(0);
        if (which & 0x1) stat += test_Gauss_Legendre_quadrature(echo);
        if (which & 0x2) stat += test_energy_mesh(echo);
        if (which & 0x4) stat += test_integrator(echo);
        return stat;
    } // all_tests

#endif // NO_UNIT_TESTS

} // namespace energy_contour
