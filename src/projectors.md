== Bloechl branch ==

get partial waves pseudized at rcut
get preliminary projectors, strictly localized inside rcut
projectors are represented in SHO basis
  - sigma is optimized
  - numax is a parameter
establish duality with SHO reconstructed projectors
  - by Gram-Schmidt orthogonalization to maintain energy ordering
use orthogonalized SHO coefficients to construct radial projector functions
create smooth partial waves with projectors as right hand side,
  - energy ordering has no ambiguity
establish duality by Gram-Schmidt again
  - overlap should already be close to unity
  - also modify the projector coefficients
  - projector coefficients are not necessarily orthogonal!

    non-local potential |~p_i> H_ij <~p_j|,      i,j < n_partial_waves[ell]
    expand in SHO functions:
    = |SHO_k> <SHO_k|~p_i> H_ij <~p_j|SHO_l> <SHO_l|,     k,l < nn_max[ell]
    n_partial_waves[ell] <= nn_max[ell]
    = |SHO_k> H_kl <SHO_l|

    density_matrix_ij = sum_n w_n <~Psi_n|~p_i> <~p_j|~Psi_n>,   Psi_n eigenstate of H
    insert SHO basis:
    density_matrix_ij = sum_n w_n <~Psi_n|SHO_k> <SHO_k|~p_i> <~p_j|SHO_l> <SHO_l|~Psi_n>
    density matrix in SHO space
    density_matrix_kl = sum_n w_n <~Psi_n|SHO_k> <SHO_l|~Psi_n>,  uses sho_projection coefficients
    transform to
    density_matrix_ij = <SHO_k|~p_i>^T density_matrix_kl <~p_j|SHO_l>^T (already implemented)

