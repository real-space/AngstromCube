#pragma once

// #include "inline_math.hxx" // set
// #include "energy_level.hxx" // TRU, SMT

namespace energy_contribution {

  int constexpr max_number=16;
  int constexpr TOTAL=0;
  int constexpr DIFFERENCE=2;
  int constexpr REFERENCE=3;
  int constexpr ELECTROSTATIC=4;
  int constexpr COULOMB=5;
  int constexpr HARTREE=6;
  int constexpr EXCHANGE_CORRELATION=7;
  int constexpr DOUBLE_COUNTING=8;
  int constexpr KINETIC=10;
  int constexpr EXTERNAL=11; // due to an external potential

  // 
  // To compute the total energy, we need
  //    + the 3D grid electrostatic energy 1/2 <Ves|rho_aug>
  //    + the 3D grid XC energy <Exc|rho>
  //    - the 3D grid XC double counting <Vxc|rho>  ?really?
  //    + the bandsum from eigenstates E_nk = <psi|H|psi>
  //    - the 3D grid double counting <Vtot|rho_valence>
  //    + the atomic TRU XC energy
  //    - the atomic SMT XC energy
  //    - the atomic TRU double counting energy
  //    + the atomic SMT double counting energy
  //    + the atomic TRU electrostatic energy (without self-interaction of the core)
  //    - the atomic TRU electrostatic energy
  //    + the atomic kinetic energy of core states
  //    + the atomic TRU kinetic energy of valence states
  //    - the atomic SMT kinetic energy of valence states
  //    ...
  //
  // Split according to 
  // TOTAL = KINETIC + EXCHANGE_CORRELATION + ELECTROSTATIC
  //    ELECTROSTATIC = COULOMB + HARTREE
  // DOUBLE_COUNTING (needed as intermediate)
  
  // Split according to CORE, SEMICORE, VALENCE, INSEPARABLE
  // Split according to TRU, SMT
  // Split atomic energy contribution after ell

  
  
//   class energy_contribution_t {
//     private:
//       static constexpr int _ellmax = 7;
//       
//     public:
//       energy_contribution_t() {
//           set(_resolved, 1 + _ellmax, 0.0);
//       } // default constructor
// 
//       double total_energy() const { return 0; }
//       inline int ellmax() const { return _ellmax; }
//       
//       void set_energy_contribution() {}
// 
//     private:
//       double _resolved[1 + _ellmax];
// 
//   }; // class energy_contribution_t

} // namespace energy_contribution
