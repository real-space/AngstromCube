#include <stdint.h> // int32_t, int8_t
#include "single_atom.h" // interface
#include <stdio.h> // printf
#include <math.h> // sqrt

int main(int argc, char *argv[]) {
    // constants
#ifndef NA
#define NA 12 // number of atoms
#endif

    // variables
    int32_t status = 0, stride = 0;
    int ia, ir2, iq;
    double Z_core[NA], sigma[NA], rcut[NA], ionic[NA], magn[NA], nve[NA], energy[NA];
    int const nr2 = 1 << 12;
    double memory[NA][nr2];
    double* pointers[NA];
    int32_t atom_id[NA], numax[NA], lmax_qlm[NA], lmax_vlm[NA];
    int8_t nn[NA][8];
    char const xc_key[] = "LDA";
    int32_t const na = NA;
    float const ar2 = 16;

//  printf("# %s:%d\n", __FILE__, __LINE__); // here

    live_atom_init_env_("control.sh", &status);
    printf("# live_atom_init_env = %d\n", status);

    live_atom_set_env_("called.from", "C", &status);
    printf("# live_atom_set_env = %d\n", status);

    for (ia = 0; ia < na; ++ia) {
        Z_core[ia] = ia + 1;
        pointers[ia] = memory[ia];
    } // ia

    live_atom_initialize_(&na, Z_core, atom_id, numax, sigma,
                        rcut, nn, ionic, magn, xc_key, &stride,
                        lmax_qlm, lmax_vlm, nve, &status);
    printf("# live_atom_initialize = %d\n", status);

    
    live_atom_get_core_density_(&na, pointers, &status);
    printf("# live_atom_get_core_density = %d\n", status);
    
    for (ia = 0; ia < 0; ++ia) {
        printf("# r, memory(r) for atom #%i\n", ia);
        for (ir2 = 0; ir2 < nr2; ++ir2) {
            printf("%.6f %.3e\n", sqrt(ir2/ar2), memory[ia][ir2]);
        } // ir2
        printf("\n");
    } // ia

    live_atom_get_compensation_charge_(&na, pointers, &status);
    printf("# live_atom_get_compensation_charge = %d\n", status);

    live_atom_set_potential_multipole_(&na, pointers, &status);
    printf("# live_atom_set_potential_multipole = %d\n", status);

    live_atom_set_density_matrix_(&na, pointers, &status);
    printf("# live_atom_set_density_matrix = %d\n", status);

    live_atom_get_hamiltonian_matrix_(&na, pointers, &status);
    printf("# live_atom_get_hamiltonian_matrix = %d\n", status);

    live_atom_get_projectors_(&na, sigma, numax, &status);
    printf("# live_atom_get_projectors = %d\n", status);

    live_atom_get_energy_contributions_(&na, energy, &status); // SEGFAULTs
    printf("# live_atom_get_energy_contributions = %d\n", status);

    live_atom_finalize_(&na, &status);
    printf("# live_atom_finalize = %d\n", status);
#undef NA
//  printf("# %s:%d\n", __FILE__, __LINE__); // here
} // main
