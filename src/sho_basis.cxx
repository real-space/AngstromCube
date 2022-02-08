#include <cstdio> // std::printf
#include <cstdlib> // std::atof, ::atoi
#include <vector> // std::vector<T>
#include <utility> // std::pair<T>, ::make_pair

#include "sho_basis.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "data_view.hxx" // view2D<T>
#include "control.hxx" // ::get
#include "xml_reading.hxx" // ::read_sequence
#include "unit_system.hxx" // Ang, _Ang
#include "sho_tools.hxx" // ::nn_max, ::nSHO_radial
#include "scattering_test.hxx" // ::expand_sho_projectors
#include "radial_grid.hxx" // radial_grid_t, ::create_radial_grid, ::destroy_radial_grid

#ifdef  HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml
  #include "rapidxml/rapidxml.hpp" // ::xml_document<>
  #include "rapidxml/rapidxml_utils.hpp" // ::file<>

  #include "xml_reading.hxx" // ::find_attribute
#endif // HAS_RAPIDXML


namespace sho_basis {
  // loads radial basis function that are expanded into a SHO basis
 
  status_t load(
        std::vector<view2D<double>> & basis
      , std::vector<int> & indirection
      , int const natoms // number of SHO basis centers
      , double const Z_core[] // Z_core[atoms]
      , int const echo // =0 log-level
  ) {
      auto const filename = control::get("sho_basis.file", "pseudo_basis.xml");
#ifndef HAS_RAPIDXML
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
              error("failed to parse %s", filename);
              return 2; // error
          } // try + catch
      } catch (...) {
          error("failed to open %s", filename);
          return 1; // error
      } // try + catch

      auto const main_node = doc.first_node("basis");
      if (!main_node) return -1; // error

      char const ellchar[] = "spdfgh?";

      auto rg = *radial_grid::create_radial_grid(256, 10.f, '='); // for display

      std::vector<std::vector<std::vector<double>>> all_coeffs;
      std::vector<double> all_sigma;
      std::vector<int> all_numax;
      std::vector<std::vector<int8_t>> all_ells, all_enns;
      
      // check all bases
      for (auto atomic = main_node->first_node("atomic"); atomic; atomic = atomic->next_sibling()) {
//        if (echo > 0) std::printf("# %s atomic= %p\n", __func__, (void*)atomic);
          auto const symbol =          (xml_reading::find_attribute(atomic, "symbol", "?"));
          auto const Z      = std::atof(xml_reading::find_attribute(atomic, "Z",     "-9"));
          auto const numax  = std::atoi(xml_reading::find_attribute(atomic, "numax", "-1"));
          auto const sigma  = std::atof(xml_reading::find_attribute(atomic, "sigma",  "0"));
          if (echo > 0) std::printf("# %s symbol= %s Z= %g numax= %d sigma= %g %s\n", __func__, symbol, Z, numax, sigma*Ang, _Ang);

          std::vector<std::vector<double>> coeffs;
          std::vector<int8_t> enns, ells;
          for (auto wave = atomic->first_node("wave"); wave; wave = wave->next_sibling()) {
              int const enn = std::atoi(xml_reading::find_attribute(wave, "n",  "0"));
              int const ell = std::atoi(xml_reading::find_attribute(wave, "l", "-1"));
              assert(ell >= 0);
              assert(enn > ell);
              size_t const n_expect = sho_tools::nn_max(numax, ell);
              auto const vec = xml_reading::read_sequence<double>(wave->value(), echo, n_expect);
              if (echo > 0) std::printf("# %s   %s-%d%c  (%ld of %ld elements)\n", __func__, symbol, enn, ellchar[ell], vec.size(), n_expect);
              assert(n_expect == vec.size());
              enns.push_back(enn);
              ells.push_back(ell);
              coeffs.push_back(vec);
          } // wave

          all_coeffs.push_back(coeffs);
          all_numax.push_back(numax);
          all_sigma.push_back(sigma);
          all_ells.push_back(ells);
          all_enns.push_back(enns);

//        if (echo > 0) std::printf("# %s atomic->next_sibling= %p\n", __func__, (void*)atomic->next_sibling());
      } // atomic

      
      for (int k = 0; k < all_coeffs.size(); ++k) {
          auto const & coeffs = all_coeffs[k];
          auto const numax  = all_numax[k];
          auto const sigma  = all_sigma[k];
          auto const & ells = all_ells[k];
          auto const & enns = all_enns[k];

          if (coeffs.size() > 0) { // plot
              assert(coeffs.size() == enns.size());
              assert(coeffs.size() == ells.size());
              if (echo > 7) {
                  int const nln = sho_tools::nSHO_radial(numax);
                  view2D<double> basis_funcs(nln, rg.n, 0.0);
                  scattering_test::expand_sho_projectors(basis_funcs.data(), basis_funcs.stride(), rg, sigma, numax);
                  std::printf("\n## radius ");
                  for (int i = 0; i < coeffs.size(); ++i) {
                      assert(ells[i] >= 0);
                      std::printf(" %d%c", enns[i], ellchar[ells[i]]);
                      assert(coeffs[i].size() == sho_tools::nn_max(numax, ells[i]));
                  } // i
                  std::printf(" :\n");
                  for (int ir = 0; ir < rg.n; ++ir) {
                      std::printf("%g", rg.r[ir]);
                      for (int i = 0; i < coeffs.size(); ++i) {
                          int const iln0 = sho_tools::ln_index(numax, ells[i], 0);
                          double wave{0};
                          for (int jrn = 0; jrn < coeffs[i].size(); ++jrn) {
                              wave += coeffs[i][jrn] * basis_funcs(iln0 + jrn,ir);
                          } // jrn
                          std::printf(" %g", wave);
                      } // i
                      std::printf("\n");
                  } // ir
                  std::printf("\n\n");
              } // echo
          } // plot
          
      } // k
      
      radial_grid::destroy_radial_grid(&rg);

      return 0;
#endif // HAS_RAPIDXML
  }

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_load(int const echo=5) {
      status_t stat(0);

      int const natoms = 1;
      std::vector<view2D<double>> basis;
      std::vector<int> indirection;
      std::vector<double> Z_core(natoms, 29.);
      stat += load(basis, indirection, natoms, Z_core.data(), echo);

      return stat;
  } // test_load
  
  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_load(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace sho_basis
