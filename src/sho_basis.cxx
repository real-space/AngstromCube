#include <cstdio> // std::printf
#include <cstdlib> // std::atof, ::atoi
#include <vector> // std::vector<T>

#include "sho_basis.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "data_view.hxx" // view2D<T>
#include "control.hxx" // ::get
#include "xml_reading.hxx" // ::read_sequence
#include "unit_system.hxx" // Ang, _Ang
#include "sho_tools.hxx" // ::nn_max

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

      // check all bases
      for (auto atomic = main_node->first_node("atomic"); atomic; atomic = atomic->next_sibling()) {
          auto const symbol =          (xml_reading::find_attribute(atomic, "symbol", "?"));
          auto const Z      = std::atof(xml_reading::find_attribute(atomic, "Z",     "-9"));
          auto const numax  = std::atoi(xml_reading::find_attribute(atomic, "numax", "-1"));
          auto const sigma  = std::atof(xml_reading::find_attribute(atomic, "sigma",  "0"));
          if (echo > 0) std::printf("# %s symbol= %s Z= %g numax= %d sigma= %g %s\n", __func__, symbol, Z, numax, sigma*Ang, _Ang);

          for (auto wave = atomic->first_node("wave"); wave; wave = wave->next_sibling()) {
              int const enn = std::atoi(xml_reading::find_attribute(wave, "n",  "0"));
              int const ell = std::atoi(xml_reading::find_attribute(wave, "l", "-1"));
              assert(ell >= 0);
              assert(enn > ell);
              int const n_expect = sho_tools::nn_max(numax, ell);
              auto const vec = xml_reading::read_sequence<double>(wave->value(), echo, n_expect);
              if (echo > 0) std::printf("# %s   %s-%d%c  (%ld of %d elements)\n", __func__, symbol, enn, ellchar[ell], vec.size(), n_expect);
              assert(n_expect == vec.size());
          } // wave

      } // atomic
      
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
      stat = load(basis, indirection, natoms, Z_core.data(), echo);

      return stat;
  } // test_load
  
  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_load(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace sho_basis
