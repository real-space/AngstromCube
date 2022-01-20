#pragma once

#include <cstdio> // std::printf
// #include <cstdint> // int64_t, int32_t, uint32_t, int8_t
#include <cassert> // assert
#include <cmath> // std::sqrt
#include <algorithm> // std::max
#include <vector> // std::vector<T>
#include <cstdlib> // std::atoi, ::atof

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer

#ifdef  HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml
  #include "rapidxml/rapidxml.hpp" // ::xml_document<>
  #include "rapidxml/rapidxml_utils.hpp" // ::file<>

  #include "xml_reading.hxx" // ::find_attribute, ::find_child
#endif // HAS_RAPIDXML

#include "xml_reading.hxx" // ::read_sequence

#include "unit_system.hxx" // eV, _eV, Ang, _Ang

#include "print_tools.hxx" // printf_vector

#include "control.hxx" // ::get

/*
 * Loading of PAW XML files generated for GPAW
 */

namespace gpaw_loading {
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_loading(int const echo=0) {
#ifndef HAS_RAPIDXML
      warn("Unable to test GPAW loading when compiled without -D HAS_RAPIDXML", 0);
      return STATUS_TEST_NOT_INCLUDED;
#else  // HAS_RAPIDXML
      SimpleTimer timer(__FILE__, __LINE__, __func__, echo);
      status_t stat(0);

      char const *filename = control::get("gpaw_loading.filename", "C.LDA");
      if (echo > 0) std::printf("# %s file=%s\n", __func__, filename);
      rapidxml::file<> infile(filename);

      rapidxml::xml_document<> doc;
      doc.parse<0>(infile.data());

      auto const paw_setup = doc.first_node("paw_setup");
      if (!paw_setup) return -1; // error
      
      auto const atom = xml_reading::find_child(paw_setup, "atom", echo);
      if (atom) {
          auto const symbol  = xml_reading::find_attribute(atom, "symbol", "?_", echo);
          auto const Z       = xml_reading::find_attribute(atom, "Z",      "-9", echo);
          auto const core    = xml_reading::find_attribute(atom, "core",    "0", echo);
          auto const valence = xml_reading::find_attribute(atom, "valence", "0", echo);
          if (echo > 0) std::printf("# %s:  <atom symbol=\"%s\" Z=\"%s\" core=\"%s\" valence=\"%s\"/>\n",
              filename, symbol, Z, core, valence);
      } else warn("no <atom> found in xml-file '%s'", filename);
      

      auto const xc_functional = xml_reading::find_child(paw_setup, "xc_functional", echo);
      if (xc_functional) {
          auto const type = xml_reading::find_attribute(xc_functional, "type", "?type");
          auto const name = xml_reading::find_attribute(xc_functional, "name", "?name");
          if (echo > 0) std::printf("# %s:  <xc_functional type=\"%s\" name=\"%s\"/>\n", filename, type, name);
      } else warn("no <xc_functional> found in xml-file '%s'", filename);

      auto const generator = xml_reading::find_child(paw_setup, "generator", echo);
      if (generator) {
          auto const type = xml_reading::find_attribute(generator, "type", "?type");
          auto const name = xml_reading::find_attribute(generator, "name", "?name");
          // caveat: generator->value() contains \n characters
          if (echo > 0) std::printf("# %s:  <generator type=\"%s\" name=\"%s\"> ... </generator>\n",
              filename, type, name);
      } else warn("no <generator> found in xml-file '%s'", filename);

      auto const ae_energy = xml_reading::find_child(paw_setup, "ae_energy", echo);
      if (ae_energy) {
          auto const kinetic       = xml_reading::find_attribute(ae_energy, "kinetic",       "0");
          auto const xc            = xml_reading::find_attribute(ae_energy, "xc",            "0");
          auto const electrostatic = xml_reading::find_attribute(ae_energy, "electrostatic", "0");
          auto const total         = xml_reading::find_attribute(ae_energy, "total",         "0");
          if (echo > 0) std::printf("# %s:  <ae_energy kinetic=\"%s\" xc=\"%s\" electrostatic=\"%s\" total=\"%s\"/>\n",
              filename, kinetic, xc, electrostatic, total);
      } else warn("no <ae_energy> found in xml-file '%s'", filename);

      auto const core_energy = xml_reading::find_child(paw_setup, "core_energy", echo);
      if (core_energy) {
          auto const kinetic = xml_reading::find_attribute(core_energy, "kinetic",       "0");
          if (echo > 0) std::printf("# %s:  <core_energy kinetic=\"%s\"/>\n", filename, kinetic);
      } else warn("no <core_energy> found in xml-file '%s'", filename);


      int ngrid{0};
      auto const radial_grid = xml_reading::find_child(paw_setup, "radial_grid", echo);
      if (radial_grid) {
          auto const eq     = xml_reading::find_attribute(radial_grid, "eq", "?eq");
          auto const a      = xml_reading::find_attribute(radial_grid, "a", ".4");
          auto const n      = xml_reading::find_attribute(radial_grid, "n", "0");
          auto const istart = xml_reading::find_attribute(radial_grid, "istart", "0");
          auto const iend   = xml_reading::find_attribute(radial_grid, "iend", "-1");
          auto const id     = xml_reading::find_attribute(radial_grid, "id", "?");
          if (echo > 0) std::printf("# %s:  <radial_grid eq=\"%s\" a=\"%s\" n=\"%s\" istart=\"%s\" iend=\"%s\" id=\"%s\"/>\n",
              filename, eq, a, n, istart, iend, id);
          ngrid = std::atoi(n); // convert to integer
      } else warn("no <radial_grid> found in xml-file '%s'", filename);

      char const *radial_state_quantities[20] = {"ae_partial_wave", "pseudo_partial_wave", "projector_function"};

      int nstates{0};
      int nwarn[] = {0, 0, 0};
      auto const valence_states = xml_reading::find_child(paw_setup, "valence_states", echo);
      if (valence_states) {
          if (echo > 0) std::printf("# %s:  <valence_states>\n", filename);
          for (auto state = valence_states->first_node(); state; state = state->next_sibling()) {
              auto const n  = xml_reading::find_attribute(state, "n",  "0");
              auto const l  = xml_reading::find_attribute(state, "l", "-1");
              auto const f  = xml_reading::find_attribute(state, "f",  "0");
              auto const rc = xml_reading::find_attribute(state, "rc", "0");
              auto const e  = xml_reading::find_attribute(state, "e",  "0");
              auto const id = xml_reading::find_attribute(state, "id", "?");
              if (echo > 0) std::printf("# %s:    <state n=\"%s\" l=\"%s\" f=\"%s\" rc=\"%s\" e=\"%s\" id=\"%s\"/>\n",
                  filename, n, l, f, rc, e, id);

              for (int iq = 0; iq < 3; ++iq) {
                  auto const q_name = radial_state_quantities[iq];
                  int radial_data_found{0};
                  for (auto child = paw_setup->first_node(); child; child = child->next_sibling()) {
                      if (0 == std::strcmp(q_name, child->name())) {
                          auto const state_id = xml_reading::find_attribute(child, "state", "?state");
                          if (0 == std::strcmp(state_id, id)) {
                              auto const grid = xml_reading::find_attribute(child, "grid", "?grid");
                              auto const vals = xml_reading::read_sequence<double>(child->value(), echo, ngrid);
                              if (echo > 0) std::printf("# %s:  <%s state=\"%s\" grid=\"%s\"> ...(%ld numbers)... </%s>\n",
                                  filename, q_name, state_id, grid, vals.size(), q_name);
                              nwarn[iq] += (vals.size() != ngrid);
                              ++radial_data_found;
                          } // state_id matches
                      } // found
                  } // attr
                  if (1 != radial_data_found) {
                      error("%s: radial state quantities in pawxml files must be defined exactly once!", filename);
                      ++stat;
                  } // found more or less than once
              } // iq
              ++nstates;

          } // state
          if (echo > 0) std::printf("# %s:  </valence_states>\n", filename);
      } else warn("no <valence_states> found in xml-file '%s'", filename);

      for (int iq = 0; iq < 3; ++iq) {
          if (nwarn[iq]) {
              warn("%s: %d %s deviate from the expected number of %d grid points",
                   filename, nwarn[iq], radial_state_quantities[iq], ngrid);
          } // warn
      } // iq

      auto const shape_function = xml_reading::find_child(paw_setup, "shape_function", echo);
      if (shape_function) {
          auto const type = xml_reading::find_attribute(shape_function, "type", "?type");
          auto const rc   = xml_reading::find_attribute(shape_function, "rc", "0");
          if (echo > 0) std::printf("# %s:  <shape_function type=\"%s\" rc=\"%s\"/>\n", filename, type, rc);
      } else warn("no <shape_function> found in xml-file '%s'", filename);

      auto const radial_quantities = {"ae_core_density", "pseudo_core_density",
          "ae_core_kinetic_energy_density", "pseudo_core_kinetic_energy_density", 
          "zero_potential"}; // , "pseudo_valence_density"};
      for (auto q_name : radial_quantities) {
          auto const q_node = xml_reading::find_child(paw_setup, q_name, echo);
          if (q_node) {
              auto const grid = xml_reading::find_attribute(q_node, "grid", "?grid");
              auto const vals = xml_reading::read_sequence<double>(q_node->value(), echo, ngrid);
              if (echo > 0) std::printf("# %s:  <%s grid=\"%s\"> ...(%ld numbers)... </%s>\n",
                  filename, q_name, grid, vals.size(), q_name);
              if (vals.size() != ngrid) warn("%s: %s has %ld but expected %d grid points",
                   filename, q_name, vals.size(), ngrid);
          } else warn("no <%s> found in xml-file '%s'", q_name, filename);
      } // q_name

      auto const kinetic_energy_differences = xml_reading::find_child(paw_setup, "kinetic_energy_differences", echo);
      if (kinetic_energy_differences) {
          auto const vals = xml_reading::read_sequence<double>(kinetic_energy_differences->value(), echo, nstates*nstates);
          if (echo > 0) std::printf("# %s:  <kinetic_energy_differences> ...(%ld numbers)... </kinetic_energy_differences>\n", filename, vals.size());
          if (vals.size() != nstates*nstates) warn("%s: expected %d^2 numbers in kinetic_energy_differences matrix but found %ld", filename, nstates, vals.size());
      } else warn("no <kinetic_energy_differences> found in xml-file '%s'", filename);

// not covert:      
//   <exact_exchange_X_matrix>  </exact_exchange_X_matrix>
//   <exact_exchange core-core="-3.462071"/>

      return stat;
#endif // HAS_RAPIDXML
  } // test_loading

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_loading(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace gpaw_loading
