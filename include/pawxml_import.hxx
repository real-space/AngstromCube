#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio>     // std::printf, ::snprintf, ::fopen, ::fprintf, ::fclose
#include <cstdint>    // int8_t
#include <cassert>    // assert
#include <cstring>    // std::strcmp
#include <cmath>      // std::sqrt
#include <algorithm>  // std::max
#include <vector>     // std::vector<T>
#include <cstdlib>    // std::atoi, ::atof

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "simple_timer.hxx" // SimpleTimer
#include "radial_grid.hxx" // ::equation_reciprocal, ::equation_exponential

#ifdef    HAS_RAPIDXML
  // git clone https://github.com/dwd/rapidxml
  #include "rapidxml.hpp" // ::xml_document<>
  #include "rapidxml_utils.hpp" // ::file<>

  #include "xml_reading.hxx" // ::find_attribute, ::find_child, ::read_sequence
#endif // HAS_RAPIDXML

#include "unit_system.hxx" // eV, _eV, Ang, _Ang
#include "print_tools.hxx" // printf_vector
#include "control.hxx" // ::get

/*
 * Loading of PAW XML files generated for GPAW
 */

namespace pawxml_import {

  struct pawxmlstate_t {
      std::vector<double> tsp[3]; // 0:true,1:smooth,2:projector function
      double e;    // energy eigenvalue
      float f, rc; // occupation, cutoff radius
      int8_t n, l; // principal quantum number, angular momentum quantum number
  }; // struct pawxmlstate_t

  struct pawxml_t {
      double Z, core, semi, valence;
      double ae_energy[4];
      double core_energy_kinetic;
      double radial_grid_a;
      double radial_grid_d;
      char   radial_grid_eq;
      int    n; // radial_grid.n
      double shape_function_rc;
      std::vector<pawxmlstate_t> states;
      std::vector<double> func[6]; // radial functions
      std::vector<double> dkin; // kinetic_energy_differences
      char xc[16];
      char Sy[4];
      status_t parse_status;
  }; // struct pawxml_t

  char const radial_state_quantities[3][20] = {"ae_partial_wave", "pseudo_partial_wave", "projector_function"};

  char const radial_quantities[6][40] = {"ae_core_density", "pseudo_core_density",
          "ae_core_kinetic_energy_density", "pseudo_core_kinetic_energy_density", 
          "zero_potential", "pseudo_valence_density"};

  inline pawxml_t parse_pawxml(char const *const filename, int const echo=0) {
#ifndef   HAS_RAPIDXML
      warn("Unable to test GPAW loading when compiled without -D HAS_RAPIDXML", 0);
      return pawxml_t();
#else  // HAS_RAPIDXML

      pawxml_t p;
      p.parse_status = __LINE__;

      if (echo > 9) std::printf("# %s file=%s\n", __func__, filename);
      rapidxml::xml_document<> doc;
      try {
          rapidxml::file<> infile(filename);
          try {
              doc.parse<0>(infile.data());
          } catch (...) {
              warn("failed to parse \"%s\"", filename);
              p.parse_status = __LINE__; return p; // error
          } // try + catch
      } catch (...) {
          warn("failed to open \"%s\"", filename);
          p.parse_status = __LINE__; return p; // error
      } // try + catch

      auto const paw_dataset = doc.first_node("paw_dataset"); // defined by ESL: https://esl.cecam.org/data/paw-xml/
      auto const paw_setup = (nullptr == paw_dataset) ? doc.first_node("paw_setup") : paw_dataset; // used by GPAW
      if (!paw_setup) { p.parse_status = __LINE__; return p; } // error

      auto const atom = xml_reading::find_child(paw_setup, "atom", echo);
      if (atom) {
          auto const symbol  = xml_reading::find_attribute(atom, "symbol", "?_", echo);
          auto const Z       = xml_reading::find_attribute(atom, "Z",      "-9", echo);
          auto const core    = xml_reading::find_attribute(atom, "core",    "0", echo);
          auto const valence = xml_reading::find_attribute(atom, "valence", "0", echo);
          if (echo > 5) std::printf("# %s:  <atom symbol=\"%s\" Z=\"%s\" core=\"%s\" valence=\"%s\"/>\n",
              filename, symbol, Z, core, valence);
          p.Z       = std::atof(Z);
          p.core    = std::atof(core);
          p.semi    = 0; // not implemented
          p.valence = std::atof(valence);
          std::snprintf(p.Sy, 3, "%s", symbol);
      } else warn("<atom> not found in xml-file '%s'", filename);

      auto const xc_functional = xml_reading::find_child(paw_setup, "xc_functional", echo);
      if (xc_functional) {
          auto const type = xml_reading::find_attribute(xc_functional, "type", "?type");
          auto const name = xml_reading::find_attribute(xc_functional, "name", "?name");
          if (echo > 5) std::printf("# %s:  <xc_functional type=\"%s\" name=\"%s\"/>\n", filename, type, name);
          std::snprintf(p.xc, 15, "%s", name);
      } else warn("<xc_functional> not found in xml-file '%s'", filename);

      auto const generator = xml_reading::find_child(paw_setup, "generator", echo);
      if (generator) {
          auto const type = xml_reading::find_attribute(generator, "type", "?type");
          auto const name = xml_reading::find_attribute(generator, "name", "?name");
          // caveat: generator->value() contains \n characters
          if (echo > 5) std::printf("# %s:  <generator type=\"%s\" name=\"%s\"> ... </generator>\n",
              filename, type, name);
      } else warn("no <generator> found in xml-file '%s'", filename);

      auto const ae_energy = xml_reading::find_child(paw_setup, "ae_energy", echo);
      if (ae_energy) {
          auto const kinetic       = xml_reading::find_attribute(ae_energy, "kinetic",       "0");
          auto const xc            = xml_reading::find_attribute(ae_energy, "xc",            "0");
          auto const electrostatic = xml_reading::find_attribute(ae_energy, "electrostatic", "0");
          auto const total         = xml_reading::find_attribute(ae_energy, "total",         "0");
          if (echo > 5) std::printf("# %s:  <ae_energy kinetic=\"%s\" xc=\"%s\" electrostatic=\"%s\" total=\"%s\"/>\n",
              filename, kinetic, xc, electrostatic, total);
          p.ae_energy[0] = std::atof(kinetic);
          p.ae_energy[1] = std::atof(xc);
          p.ae_energy[2] = std::atof(electrostatic);
          p.ae_energy[3] = std::atof(total);
      } else warn("<ae_energy> not found in xml-file '%s'", filename);

      auto const core_energy = xml_reading::find_child(paw_setup, "core_energy", echo);
      if (core_energy) {
          auto const kinetic = xml_reading::find_attribute(core_energy, "kinetic",       "0");
          if (echo > 5) std::printf("# %s:  <core_energy kinetic=\"%s\"/>\n", filename, kinetic);
          p.core_energy_kinetic = std::atof(kinetic);
      } else warn("<core_energy> not found in xml-file '%s'", filename);

      p.n = 0;
      auto const radial_grid = xml_reading::find_child(paw_setup, "radial_grid", echo);
      if (radial_grid) {
          auto const eq     = xml_reading::find_attribute(radial_grid, "eq", "?equation");
          auto const a      = xml_reading::find_attribute(radial_grid, "a", ".4");
          auto const n      = xml_reading::find_attribute(radial_grid, "n", "0");
          auto const istart = xml_reading::find_attribute(radial_grid, "istart", "0");
          auto const iend   = xml_reading::find_attribute(radial_grid, "iend", "-1");
          auto const id     = xml_reading::find_attribute(radial_grid, "id", "?");
          auto const d      = xml_reading::find_attribute(radial_grid, "d", "1");
          if (echo > 5) std::printf("# %s:  <radial_grid eq=\"%s\" a=\"%s\" n=\"%s\" istart=\"%s\" iend=\"%s\" id=\"%s\"/>\n",
              filename, eq, a, n, istart, iend, id);
          p.n = std::atoi(n);
          p.radial_grid_a = std::atof(a);
          p.radial_grid_d = std::atof(d);
          p.radial_grid_eq = *eq;
          if (std::string(eq) == "r=a*i/(n-i)") {
              p.radial_grid_eq = radial_grid::equation_reciprocal;
          } else if (std::string(eq) == "r=a*(exp(d*i)-1)") {
              p.radial_grid_eq = radial_grid::equation_exponential;
          } else {
              warn("%s: radial grid neither exponential nor reciprocal, found %s", filename, eq);
          }
          if (0 != std::atoi(istart))                error("%s: assume a radial grid starting at 0", filename);
          if (p.n != std::atoi(iend) + 1)            error("%s: assume a radial grid starting from 0 to n-1", filename);
      } else warn("<radial_grid> not found in xml-file '%s'", filename);

      int istate{0};
      p.states.resize(0);
      int nwarn[] = {0, 0, 0};
      auto const valence_states = xml_reading::find_child(paw_setup, "valence_states", echo);
      if (valence_states) {
          if (echo > 5) std::printf("# %s:  <valence_states>\n", filename);
          for (auto state = valence_states->first_node(); state; state = state->next_sibling()) {
              auto const n  = xml_reading::find_attribute(state, "n",  "0");
              auto const l  = xml_reading::find_attribute(state, "l", "-1");
              auto const f  = xml_reading::find_attribute(state, "f",  "0");
              auto const rc = xml_reading::find_attribute(state, "rc", "0");
              auto const e  = xml_reading::find_attribute(state, "e",  "0");
              auto const id = xml_reading::find_attribute(state, "id", "?");
              if (echo > 5) std::printf("# %s:    <state n=\"%s\" l=\"%s\" f=\"%s\" rc=\"%s\" e=\"%s\" id=\"%s\"/>\n",
                  filename, n, l, f, rc, e, id);

              pawxmlstate_t s;
              s.n  = std::atoi(n);
              s.l  = std::atoi(l);
              s.f  = std::atof(f);
              s.rc = std::atof(rc);
              s.e  = std::atof(e);

              for (int iq = 0; iq < 3; ++iq) {
                  auto const q_name = radial_state_quantities[iq];
                  int radial_data_found{0};
                  for (auto child = paw_setup->first_node(); child; child = child->next_sibling()) {
                      if (0 == std::strcmp(q_name, child->name())) {
                          auto const state_id = xml_reading::find_attribute(child, "state", "?state");
                          if (0 == std::strcmp(state_id, id)) {
                              auto const grid = xml_reading::find_attribute(child, "grid", "?grid");
                              auto const vals = xml_reading::read_sequence<double>(child->value(), echo, p.n);
                              if (echo > 8) std::printf("# %s:  <%s state=\"%s\" grid=\"%s\"> ...(%ld numbers)... </%s>\n",
                                  filename, q_name, state_id, grid, vals.size(), q_name);
                              nwarn[iq] += (vals.size() != p.n);
                              s.tsp[iq] = vals;
                              ++radial_data_found;
                          } // state_id matches
                      } // found
                  } // attr
                  if (1 != radial_data_found) error("%s: radial state quantities in pawxml files must be defined exactly once!", filename);
              } // iq
              p.states.push_back(s);
              ++istate;

          } // state
          if (echo > 5) std::printf("# %s:  </valence_states>\n", filename);
      } else warn("<valence_states> not found in xml-file '%s'", filename);
      int const nstates = istate;
      assert(nstates == p.states.size());

      for (int iq = 0; iq < 3; ++iq) {
          if (nwarn[iq]) {
              warn("%s: %d %s deviate from the expected number of %d grid points",
                   filename, nwarn[iq], radial_state_quantities[iq], p.n);
          } // nwarn
      } // iq

      auto const shape_function = xml_reading::find_child(paw_setup, "shape_function", echo);
      if (shape_function) {
          auto const type = xml_reading::find_attribute(shape_function, "type", "?type");
          auto const rc   = xml_reading::find_attribute(shape_function, "rc", "0");
          if (echo > 5) std::printf("# %s:  <shape_function type=\"%s\" rc=\"%s\"/>\n", filename, type, rc);
          if (std::strcmp(type, "gauss")) error("%s: assume a shape_function type gauss", filename);
          p.shape_function_rc = std::atof(rc);
      } else warn("<shape_function> not found in xml-file '%s'", filename);

      for (int iq = 0; iq < 5; ++iq) {
          auto const q_name = radial_quantities[iq];
          auto const q_node = xml_reading::find_child(paw_setup, q_name, echo);
          if (q_node) {
              auto const grid = xml_reading::find_attribute(q_node, "grid", "?grid");
              auto const vals = xml_reading::read_sequence<double>(q_node->value(), echo, p.n);
              if (echo > 7) std::printf("# %s:  <%s grid=\"%s\"> ...(%ld numbers)... </%s>\n",
                  filename, q_name, grid, vals.size(), q_name);
              if (vals.size() != p.n) warn("%s: %s has %ld but expected %d grid points",
                   filename, q_name, vals.size(), p.n);
              p.func[iq] = vals;
          } else warn("<%s> not found in xml-file '%s'", q_name, filename);
      } // iq

      auto const kinetic_energy_differences = xml_reading::find_child(paw_setup, "kinetic_energy_differences", echo);
      if (kinetic_energy_differences) {
          auto const vals = xml_reading::read_sequence<double>(kinetic_energy_differences->value(), echo, nstates*nstates);
          if (echo > 7) std::printf("# %s:  <kinetic_energy_differences> ...(%ld numbers)... </kinetic_energy_differences>\n", filename, vals.size());
          if (vals.size() != nstates*nstates) warn("%s: expected %d^2 numbers in kinetic_energy_differences matrix but found %ld", filename, nstates, vals.size());
          p.dkin = vals;
      } else warn("<kinetic_energy_differences> not found in xml-file '%s'", filename);

// not covered:
//   <exact_exchange_X_matrix>  </exact_exchange_X_matrix>
//   <exact_exchange core-core="-3.462071"/>

      p.parse_status = 0;
      return p;
#endif // HAS_RAPIDXML
  } // parse_pawxml



#ifdef    NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  inline status_t test_loading(int const echo=0) {
      SimpleTimer timer(__FILE__, __LINE__, __func__, echo);

      auto const filename = control::get("pawxml_import.test.filename", "C.LDA");
      if (echo > 5) std::printf("# pawxml_import.test.filename=%s\n", filename);
      auto const p = parse_pawxml(filename, echo);

      char const ellchar[8] = "spdfghi";
      int const repeat_file = control::get("pawxml_import.test.repeat", 0.);
      if (echo > 5) std::printf("# pawxml_import.test.repeat=%d\n", repeat_file);
      if (repeat_file) {
          // create an XML file that should reproduce this
          char outfile[992]; std::snprintf(outfile, 992, "%s.repeat.xml", filename);
          auto const f = std::fopen(outfile, "w");
          std::fprintf(f, "<?xml version=\"1.0\"?>\n"
              "<paw_setup version=\"0.6\">\n"
              "\t<!-- Element setup for the Projector Augmented Wave method -->\n"
              "\t<!-- Units: Hartree and Bohr radii.                        -->\n");
          std::fprintf(f, "\t<atom symbol=\"%s\" Z=\"%g\" core=\"%g\" valence=\"%g\"/>\n", p.Sy, p.Z, p.core, p.valence);
          std::fprintf(f, "\t<xc_functional name=\"%s\"/>\n", p.xc);
          std::fprintf(f, "\t<ae_energy kinetic=\"%.6f\" xc=\"%.6f\" electrostatic=\"%.6f\" total=\"%.6f\"/>\n",
                              p.ae_energy[0], p.ae_energy[1], p.ae_energy[2], p.ae_energy[3]);
          std::fprintf(f, "\t<core_energy kinetic=\"%.6f\"/>\n", p.core_energy_kinetic);
          std::fprintf(f, "\t<valence_states>\n");
          for (auto & s : p.states) {
              char id[8]; std::snprintf(id, 8, "%s-%d%c", p.Sy, s.n, ellchar[s.l]);
              std::fprintf(f, "\t\t<state n=\"%d\" l=\"%d\" f=\"%g\" rc=\"%g\" e=\"%g\" id=\"%s\"/>\n",
                                        s.n,     s.l,     s.f,     s.rc,     s.e,       id);
          } // states
          std::fprintf(f, "\t</valence_states>\n");
          std::fprintf(f, "\t<radial_grid eq=\"r=a*i/(n-i)\" a=\"%g\" n=\"%d\" istart=\"0\" iend=\"%d\" id=\"g1\"/>\n",
                                               p.radial_grid_a,     p.n,                     p.n - 1);
          std::fprintf(f, "\t<shape_function type=\"gauss\" rc=\"%.12e\"/>\n", p.shape_function_rc);

          // radial_quantities
          for (int iq = 0; iq < 5; ++iq) {
              std::fprintf(f, "\t<%s grid=\"g1\">\n", radial_quantities[iq]);
              for (auto val : p.func[iq]) {
                  std::fprintf(f, " %g", val);
              } // val
              std::fprintf(f, "\n\t</%s>\n", radial_quantities[iq]);
          } // iq

          // radial_state_quantities
          for (auto const & s : p.states) {
              char id[8]; std::snprintf(id, 7, "%s-%d%c", p.Sy, s.n, ellchar[s.l]);
              for (int iq = 0; iq < 3; ++iq) {
                  std::fprintf(f, "\t<%s state=\"%s\" grid=\"g1\">\n", radial_state_quantities[iq], id);
                  for (auto val : s.tsp[iq]) {
                      std::fprintf(f, " %g", val);
                  } // val
                  std::fprintf(f, "\n\t</%s>\n", radial_state_quantities[iq]);
              } // iq
          } // s

          { // kinetic energy differences
              std::fprintf(f, "\t<kinetic_energy_differences>");
              int const n = p.states.size();
              for (int i = 0; i < n; ++i) {
                  for (int j = 0; j < n; ++j) {
                      std::fprintf(f, "%s%g", j?" ":"\n\t\t", p.dkin[i*n + j]);
                  } // j
              } // i
              std::fprintf(f, "\n\t</kinetic_energy_differences>\n");
          } // scope

          std::fprintf(f, "</paw_setup>\n");
          std::fclose(f);
      } // repeat_file

      int const plot_file = control::get("pawxml_import.test.plot", double(echo > 7));
      if (echo > 5) std::printf("# pawxml_import.test.plot=%d   # 0:none 1:true"
                                " 2:smooth 4:projectors 7:all\n", plot_file);
      int const n = p.states.size();
      char const what[3][32] = {"true partial waves", "smooth partial waves", "projectors"};
      for (int iq = 0; iq < 3; ++iq) { // TRU=0, SMT=1, PRJ=2
          if (plot_file & (1 << iq)) { // query bits
              std::printf("\n\n# plot %s from %s in xmgrace -nxy\n", what[iq], filename);
              std::printf("# radius/Bohr");
              for (int i = 0; i < n; ++i) {
                  std::printf(", %d%c-projector", p.states[i].n, ellchar[p.states[i].l]);
              } // i
              std::printf("\n");
              for (int ir = 1; ir < p.n; ++ir) {
                  double const r = p.radial_grid_a*
                    ((radial_grid::equation_exponential == p.radial_grid_eq) ? (std::exp(p.radial_grid_d*ir) - 1) : (ir/double(p.n - ir)));
                  std::printf("%.6f", r);
                  for (int i = 0; i < n; ++i) {
                      std::printf(" %g", p.states[i].tsp[iq][ir]); // the state's radial function
                  } // i
                  std::printf("\n");
              } // ir
              std::printf("\n\n");
          } // bit is set
      } // iq

      return p.parse_status;
  } // test_loading

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_loading(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace pawxml_import
