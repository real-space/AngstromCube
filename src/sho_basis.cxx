// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdlib> // std::atof, ::atoi
#include <cstdint> // uint8_t
#include <cassert> // assert
#include <vector> // std::vector<T>
#include <map> // std::map<Key,T>

#include "sho_basis.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "data_view.hxx" // view2D<T>
#include "control.hxx" // ::get
#include "xml_reading.hxx" // ::read_sequence
#include "unit_system.hxx" // Ang, _Ang
#include "sho_tools.hxx" // ::nn_max, ::nSHO_radial, ::nSHO, ::order_*
#include "scattering_test.hxx" // ::expand_sho_projectors
#include "radial_grid.hxx" // ::create_radial_grid, ::destroy_radial_grid, ::equation_equidistant
#include "complex_tools.hxx" // complex_name
#include "sho_unitary.hxx" // ::Unitary_SHO_Transform
#include "print_tools.hxx" // printf_vector

#ifdef    HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml
  #include "rapidxml.hpp" // ::xml_document<>
  #include "rapidxml_utils.hpp" // ::file<>

  #include "xml_reading.hxx" // ::find_attribute
#endif // HAS_RAPIDXML


namespace sho_basis {
  // loads radial basis function that are expanded into a SHO basis

  struct RadialFunction {
     std::vector<double> vec;
     int8_t enn = -1, ell = -1;
  }; // RadialFunction

  struct RadialFunctionSet {
     std::vector<RadialFunction> vec;
     double sigma = 0;
     int numax = -1;
  }; // RadialFunctionSet

  struct SpeciesSet {
     std::map<uint8_t,RadialFunctionSet> map;
     double Z_core;
     char symbol[8];
     int numax_min = 999;
     int numax_max =  -1;
  }; // SpeciesSet


  char ellchar(int const ell) {
      char const _ellchar[8] = "spdfgh?";
      return _ellchar[ell & 0x7];
  } // ellchar

  double norm2(std::vector<double> const & vec) {
      double n2{0};
      for (auto v : vec) {
          n2 += v*v;
      } // v
      return n2;
  } // norm2

  status_t load(
        RadialFunctionSet const* & rfset // result
      , double const Z_core
      , int const numax_in // =-1
      , int const echo // =0 log-level
      , bool const plot=false
  ) {

      static std::map<double,SpeciesSet> _map; // map_Key=Z_core

      rfset = nullptr;
  if (0 == _map.count(Z_core)) { // scope: load file into static _map

      // load basis data from file
      auto const filename = control::get("sho_basis.file", "pseudo_basis.xml");
#ifndef   HAS_RAPIDXML
      warn("Unable load sho_basis from \'%s\' when compiled without -D HAS_RAPIDXML", filename);
      return STATUS_TEST_NOT_INCLUDED;
#else  // HAS_RAPIDXML

      if (echo > 9) std::printf("# %s file=%s\n", __func__, filename);
      rapidxml::xml_document<> doc;
      try {
          rapidxml::file<> infile(filename);
          try {
              doc.parse<0>(infile.data());
          } catch (...) {
              warn("failed to parse \"%s\"", filename);
              return 2; // error
          } // try + catch
      } catch (...) {
          warn("failed to open \"%s\"", filename);
          return 1; // error
      } // try + catch

      auto const main_node = doc.first_node("basis");
      if (!main_node) return -1; // error

      double max_norm_err{0};

      // check all bases
      int n_species{0};
      for (auto species = main_node->first_node("species"); species; species = species->next_sibling()) {
//        if (echo > 0) std::printf("# %s species= %p\n", __func__, (void*)species);
          auto const symbol =      xml_reading::find_attribute(species, "symbol", "?");
          auto const Z = std::atof(xml_reading::find_attribute(species, "Z", "-9"));

          auto & ss = _map[Z]; ss.Z_core = Z; std::snprintf(ss.symbol, 8, "%s", symbol);

          int n_sets{0};
          for (auto set = species->first_node("set"); set; set = set->next_sibling()) {
              auto const numax = std::atoi(xml_reading::find_attribute(set, "numax", "-1"));
              auto const sigma = std::atof(xml_reading::find_attribute(set, "sigma",  "0"));
              if (echo > 5) std::printf("\n# %s symbol= %s Z= %g numax= %d sigma= %g %s\n", __func__, symbol, Z, numax, sigma*Ang, _Ang);

              auto & rfs = ss.map[numax];
              rfs.numax = numax;
              rfs.sigma = sigma;
              ss.numax_min = std::min(ss.numax_min, numax);
              ss.numax_max = std::max(ss.numax_max, numax);

              int n_waves{0}, n_basis{0};
              for (auto wave = set->first_node("wave"); wave; wave = wave->next_sibling()) {
                  auto const enn = std::atoi(xml_reading::find_attribute(wave, "n",  "0"));
                  auto const ell = std::atoi(xml_reading::find_attribute(wave, "l", "-1"));
                  assert(ell >= 0);
                  assert(enn > ell);
                  RadialFunction rf;
                  rf.ell = ell;
                  rf.enn = enn;
                  size_t const n_expect = sho_tools::nn_max(numax, ell);
                  rf.vec = xml_reading::read_sequence(wave->value(), echo, n_expect);
                  auto const norm_err = norm2(rf.vec) - 1;
                  max_norm_err = std::max(max_norm_err, std::abs(norm_err));
                  if (echo > 7) std::printf("# %s   %s-%d%c  (%ld of %ld elements), normalization %.1e\n", __func__, symbol, enn, ellchar(ell), rf.vec.size(), n_expect, norm_err);
                  assert(n_expect == rf.vec.size());
                  rfs.vec.push_back(rf);
                  ++n_waves;
                  n_basis += (2*ell + 1);
              } // wave
              if (echo > 3) std::printf("# %s symbol= %s Z= %g numax= %d sigma= %g %s has %d radial, %d basis functions\n",
                                            __func__, symbol, Z, numax, sigma*Ang, _Ang, n_waves, n_basis);
              ++n_sets;
          } // set
          if (echo > 3) std::printf("# %s symbol= %s Z= %g contains %d sets\n", __func__, symbol, Z, n_sets);

          ++n_species;
      } // species
      if (echo > 1) std::printf("# file \'%s\' contains %d species, largest normalization error is %.1e\n", filename, n_species, max_norm_err);

      // Idea: create a coefficient matrix that holds the wave functions coefficients.
      //       Modify the module sho_hamiltonian to incorporate such a matrix
      //       which defaults to unity if no basis file is specified.
      //       Perform the solving in that basis to make a cheap method possible
      //       e.g. for molecules.
      //       The coefficient matrix includes the Cartesian to Radial sho_transform so we can label the 
      //       basis functions with sharp n,l,m-quantum numbers


  } // scope: load file into static _map

      auto const found_Z = _map.count(Z_core); // try again
      if (1 == found_Z) {
          auto & speciesset = _map[Z_core];
          if (echo > 3) std::printf("# %s found Z= %g --> %s %g\n", __func__, Z_core, speciesset.symbol, speciesset.Z_core);
          assert(1 == found_Z && "std::map can only return 1 or 0 on count");
          uint8_t nu = 0;
          if (numax_in < 0) {
              nu = speciesset.numax_min; // use minimum basis
              if (echo > 3) std::printf("# %s found Z= %g --> minimum numax= %d\n", __func__, Z_core, nu);
          } else if (numax_in > speciesset.numax_max) {
              nu = speciesset.numax_max; // use maximum basis
              if (echo > 3) std::printf("# %s found Z= %g --> maximum numax= %d\n", __func__, Z_core, nu);
          } else {
              nu = uint8_t(numax_in); // use exactly the SHO basis size requested
          }
          auto const found_nu = speciesset.map.count(nu);
          if (found_nu) {
              rfset = & speciesset.map[nu];
              assert(1 == found_nu && "std::map can only return 1 or 0 on count");
              if (echo > 3) std::printf("# %s found Z= %g numax= %d ptr= %p\n", __func__, Z_core, nu, (void*)rfset);
  
              if (plot) { // scope: plot the basis functions of the requsted species
                  std::printf("# %s plot for Z= %g numax= %d\n", __func__, Z_core, nu, (void*)rfset);
                  auto rg = *radial_grid::create_radial_grid(256, 10.f, radial_grid::equation_equidistant); // for display
                  int const numax = rfset->numax;
                  int const nln = sho_tools::nSHO_radial(numax);
                  view2D<double> basis_funcs(nln, rg.n, 0.);
                  scattering_test::expand_sho_projectors(basis_funcs.data(), basis_funcs.stride(), rg, rfset->sigma, numax);
                  std::printf("\n## radius ");
                  auto const & rfs = rfset->vec;
                  for (auto & rf : rfs) {
                      assert(rf.ell >= 0 && "ell unphysical");
                      std::printf(" %d%c", rf.enn, ellchar(rf.ell));
                      assert(rf.vec.size() == sho_tools::nn_max(numax, rf.ell) && "inconsistent number of coefficients");
                  } // rf
                  std::printf(" (numax= %d)\n", numax);
                  for (int ir = 0; ir < rg.n; ++ir) {
                      std::printf("%g ", rg.r[ir]);
                      for (auto & rf : rfs) {
                          int const iln0 = sho_tools::ln_index(numax, rf.ell, 0);
                          double wave{0};
                          for (int jrn = 0; jrn < rf.vec.size(); ++jrn) {
                              wave += rf.vec[jrn] * basis_funcs(iln0 + jrn,ir);
                          } // jrn
                          std::printf(" %g", wave);
                      } // rf
                      std::printf("\n");
                  } // ir
                  std::printf("\n\n");
                  radial_grid::destroy_radial_grid(&rg);
              } // plot

              return 0;
          } else {
              warn("numax= %d for species Z= %g not found", nu, Z_core);
              assert(0 == found_nu && "std::map can only return 1 or 0 on count");
              return 3;
          }
      } else {  
          warn("species Z= %g not found", Z_core);
          assert(0 == found_Z && "std::map can only return 1 or 0 on count");
          return 4;
      } // found_Z

#endif // HAS_RAPIDXML
  } // load



  status_t get(double & sigma, int & numax, int & nbasis, double const Z_core, int const echo) {
      RadialFunctionSet const* rfset{nullptr};
      auto const numax_in = numax;
      auto const load_stat = load(rfset, Z_core, numax, echo);
      if (0 == load_stat) {
          assert(nullptr != rfset);
          sigma = rfset->sigma;
          numax = rfset->numax;
          assert(numax >= 0);
          auto const & rfs = rfset->vec;
          nbasis = 0;
          for (auto rf : rfs) {
             nbasis += (2*rf.ell + 1);
          } // rf

          if (echo > 3) {
              std::printf("# found basis set for Z= %g numax= %d --> %d sigma= %g Bohr has %ld radial, %d basis functions\n",
                                       Z_core, numax_in, numax, sigma, rfs.size(), nbasis);
              std::printf("# matrix for Z= %g will be %d x %d\n", Z_core, sho_tools::nSHO(numax), nbasis);
          } // echo
      } else {
          if (echo > 1) std::printf("# loading for Z= %g numax= %d failed, load status= %i\n",
                                       Z_core, numax_in, int(load_stat));
      } // loading successful
      return load_stat;
  } // get


  template <typename complex_t>
  status_t generate(
        view2D<complex_t> & matrix // result: on successful exit this is a nSHO(numax) x nbasis matrix
      , double & sigma
      , int & numax // SHO basis size parameter
      , double const Z_core // nuclear charge
      , int const echo // =0, verbosity
  ) {
      if (echo > 1) std::printf("# %s<%s>(Z= %g, numax= %d)\n", __func__, complex_name<complex_t>(), Z_core, numax);

      // make sure that the pseudo_basis.xml file has been loaded for this Z
      RadialFunctionSet const * rfset = nullptr;
      auto const load_stat = load(rfset, Z_core, numax, echo); // should return a RadialFunctionSet
      if (echo > 3) std::printf("# loading status= %i ptr= %p\n", int(load_stat), (void*)rfset);
      // get nbasis
      if (0 != load_stat) {
          warn("failed to load a pseudo_basis for Z= %g numax= %d", Z_core, numax);
          return load_stat;
      } // loading failed
      if (nullptr == rfset) {
          warn("loading a pseudo_basis for Z= %g numax= %d produced nullptr", Z_core, numax);
          return 1;
      } else {
          sigma = rfset->sigma;
          numax = rfset->numax;
          auto const & rfs = rfset->vec;
          if (echo > 5) std::printf("# found for Z= %g numax= %d sigma= %g Bohr\n", Z_core, numax, sigma);

          assert(numax >= 0);
          auto const nsho = sho_tools::nSHO(numax);

          int nbasis{0};
          for (auto rf : rfs) {
             nbasis += (2*rf.ell + 1);
          } // rf
          if (echo > 7) std::printf("# basis for Z= %g numax= %d has %ld radial, %d basis functions\n", Z_core, numax, rfs.size(), nbasis);
          if (echo > 5) std::printf("# basis matrix for Z= %g numax= %d will be %d x %d\n", Z_core, numax, nsho, nbasis);

          if (nbasis > nsho) {
              error("there are more basis functions than SHO functions, %d > %d", nbasis, nsho);
              return -1;
          } // basis is larger than number of SHO functions

          sho_unitary::Unitary_SHO_Transform sho_transform(numax, echo);

          view2D<double> sho_matrix(nsho, nsho, 0.0);
          sho_transform.construct_dense_matrix(sho_matrix.data(), numax, sho_matrix.stride(), sho_tools::order_zyx, sho_tools::order_Elnm);

          if (echo > 9) {
              std::printf("# sho_matrix (%d x %d):\n", nsho, nsho);
              for (int i = 0; i < nsho; ++i) {
                  std::printf("#%3i  ", i);
                  for (int j = 0; j < nsho; ++j) {
                      std::printf(" %g", sho_matrix(i,j));
                  } // j
                  std::printf("\n");
              } // i
              std::printf("\n");
          } // echo

          std::vector<uint8_t> enn_basis(nbasis, 1), ell_basis(nbasis, 0), rf_index(nbasis, 0);
          std::vector<int8_t> emm_basis(nbasis);
          { // scope: fill quantum number arrays

              if (echo > 5) std::printf("# Z= %g has", Z_core);
              int ibasis{0}, iradial{0};
              for (auto rf : rfs) {
                  if (1) {
                      auto const n_expect = sho_tools::nn_max(numax, rf.ell);
                      assert(n_expect == rf.vec.size());
                  } // sanity checks
                  for (int emm = -rf.ell; emm <= rf.ell; ++emm) {
                        enn_basis[ibasis] = rf.enn;
                        ell_basis[ibasis] = rf.ell;
                        emm_basis[ibasis] = emm;
                        rf_index[ibasis] = iradial;
                        ++ibasis;
                  } // emm
                  if (echo > 5) std::printf(" %d%c", rf.enn, ellchar(rf.ell)); 
                  ++iradial;
              } // rf
              assert(nbasis == ibasis);
              assert(rfs.size() == iradial);
              if (echo > 5) std::printf(", %ld radial functions\n", rfs.size()); 
              
              if (echo > 9) {
                  std::printf("# enn_basis "); printf_vector(" %2d", enn_basis);
                  std::printf("# ell_basis "); printf_vector(" %2d", ell_basis);
                  std::printf("# emm_basis "); printf_vector(" %2d", emm_basis);
                  std::printf("# rf_index  "); printf_vector(" %2d", rf_index);
              } // echo
          } // scope

          // combine sho_matrix and radial function coefficients
          matrix = view2D<complex_t>(nsho, nbasis, complex_t(0)); // get memory and initialize zero
          for (int isho = 0; isho < nsho; ++isho) { // SHO basis functions in order_zyx
              for (int j = 0; j < nbasis; ++j) {
                  double c{0};
                  int const nn = sho_tools::nn_max(numax, ell_basis[j]);
                  auto const & coeff = rfs[rf_index[j]].vec;
                  for (int krn = 0; krn < nn; ++krn) { // contract over radial SHO basis functions
                      int const ksho = sho_tools::Elnm_index(ell_basis[j], krn, emm_basis[j]);
                      c += sho_matrix(isho,ksho)*coeff[krn];
                  } // krn
                  matrix(isho,j) = c;
              } // j (basis)
          } // isho (order_zyx)

          if (echo > 9) {
              std::printf("# matrix (%d x %d):\n", nsho, nbasis);
              for (int isho = 0; isho < nsho; ++isho) {
                  std::printf("#%3i  ", isho);
                  for (int j = 0; j < nbasis; ++j) {
                      std::printf(" %g", std::real(matrix(isho,j))); // these numbers will never have an imaginary part
                  } // j
                  std::printf("\n");
              } // i
              std::printf("\n");
          } // echo

      } // rfset

      if (echo > 3) std::printf("# %s done\n\n", __func__);

      return load_stat;
  } // generate

  // explicit template instantiations
  template status_t generate<std::complex<double>>(view2D<std::complex<double>> &, double &, int &, double, int);
  template status_t generate<std::complex<float >>(view2D<std::complex<float >> &, double &, int &, double, int);
  template status_t generate<double              >(view2D<double              > &, double &, int &, double, int);
  template status_t generate<float               >(view2D<float               > &, double &, int &, double, int);


#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_load(int const echo=5) {
      auto const Z = control::get("sho_basis.test.Z", 29.);
      status_t stat(0);
      RadialFunctionSet const* rfset;
      int const load_echo = control::get("sho_basis.load.echo", 0.);
      stat += load(rfset, Z, -1, load_echo); // load for the first time
      int const nu_min = control::get("sho_basis.plot.min", 5.);
      int const nu_max = control::get("sho_basis.plot.max", 4.);
      for (int numax = nu_min; numax <= nu_max; ++numax) {
          bool constexpr plot = true;
          load(rfset, Z, numax, echo, plot); // ignore status for non-existing numax entries
      } // numax
      return stat;
  } // test_load

  status_t test_generate(int const echo=5) {
      double const Z = control::get("sho_basis.test.Z", 29.);
      int numax      = control::get("sho_basis.test.numax", -1.);
      status_t stat(0);
      double sigma;
      { view2D<std::complex<double>> m; stat += generate(m, sigma, numax, Z, echo  ); }
      { view2D<std::complex<float >> m; stat += generate(m, sigma, numax, Z, echo/8); }
      { view2D<float               > m; stat += generate(m, sigma, numax, Z, echo/8); }
      { view2D<double              > m; stat += generate(m, sigma, numax, Z, echo/8); }
      return stat;
  } // test_generate

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_load(echo);
      stat += test_generate(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace sho_basis
