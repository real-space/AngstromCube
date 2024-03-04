#! /usr/bin/env julia
######################################################
### Julia interface for the LiveAtom library
######################################################
### when running for the first time you may need to activate the following line
# import Pkg; Pkg.add("Formatting")

const LA = "./libliveatom.so"

plot = false

status = zeros(Int32, 1) # int32_t *status is the last argument of all calls to LiveAtom

@ccall LA.live_atom_init_env_("control.sh"::Cstring, status::Ref{Int32})::Cvoid
# ccall(("live_atom_init_env_", LA), Cvoid, (Cstring, Ref{Int32}), "control.sh", status)
println("### live_atom_init_env_ = ", status[1])

@ccall LA.live_atom_set_env_("called.from"::Cstring, "Julia"::Cstring, status::Ref{Int32})::Cvoid
# ccall(("live_atom_set_env_", LA), Cvoid, (Cstring, Cstring, Ref{Int32}), "called.from", "Julia", status)
println("### live_atom_init_env = ", status[1])

# the total number of atoms
const natoms = 10
const ar2 = 16.
const nr2 = 2^12
pointers = Vector{Vector{Float64}}(undef, natoms) # this can be passed to double** arguments

# prepare arrays for atom properties
Z_core = ones(Float64, natoms)
atom_id = zeros(Int32, natoms)
numax = ones(Int32, natoms)
sigma = ones(Float64, natoms)
rcut = ones(Float64, natoms)
nn = zeros(Int8, 8, natoms)
ionization = zeros(Float64, natoms)
magnetization = zeros(Float64, natoms)
xc_key = "LDA"
stride = 0
lmax_qlm = ones(Int32, natoms)
lmax_vlm = ones(Int32, natoms)
n_valence_electrons = zeros(Float64, natoms)

# numax .*= -9 # if we want to signal to load pawxml files

for ia = 1:natoms
    Z_core[ia] = ia
    atom_id[ia] = ia - 1
    pointers[ia] = zeros(Float64, nr2)
end # ia
# initialize the atoms using pot/Zeff.00Z files
@ccall LA.live_atom_initialize_(
              natoms::Ref{Int32}
            , Z_core::Ref{Float64}
            , atom_id::Ref{Int32}
            , numax::Ref{Int32}
            , sigma::Ref{Float64}
            , rcut::Ref{Float64}
            , nn::Ref{Int8}
            , ionization::Ref{Float64}
            , magnetization::Ref{Float64}
            , xc_key::Cstring
            , stride::Ref{Int32}
            , lmax_qlm::Ref{Int32}
            , lmax_vlm::Ref{Int32}
            , n_valence_electrons::Ref{Float64}
            , status::Ref{Int32}
            )::Cvoid
println("### live_atom_initialize_ = ", status[1])


# get the smooth core density
@ccall LA.live_atom_get_core_density_(natoms::Ref{Int32}, pointers::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### live_atom_get_core_density_ = ", status[1])
if plot
    for ia = 1:natoms
        core_density = pointers[ia]
        println("### r, core_density(r) for atom #", ia)
        for ir2 = 1:nr2
            println(sqrt((ir2 - 1)/ar2), " ", core_density[ir2])
        end # ir2
        println("")
    end # ia
end # plot

if true
    # get start waves
    waves = zeros(Float64, 4096, 6, natoms)
    occupation = zeros(Float64, 2, 6, natoms)
    @ccall LA.live_atom_get_start_waves_(natoms::Ref{Int32}, waves::Ref{Float64}, occupation::Ref{Float64}, status::Ref{Int32})::Cvoid
    println("### get_start_waves_ = ", status[1])
    if plot
        for ia = 1:natoms
            println("### r, waves(r) for atom #", ia, " occupied ", join(occupation[1,:,ia] + occupation[2,:,ia], ' '))
            for ir2 = 1:nr2
                println(sqrt((ir2 - 1)/ar2), " ", join(waves[ir2,:,ia], ' '))
            end # ir2
            println("")
        end # ia
    end # plot
end # true

# set density matrix
atom_rho = Vector{Vector{Float64}}(undef, natoms)
for ia = 1:natoms
    atom_rho[ia] = zeros(Float64, nr2)
end # ia
@ccall LA.live_atom_set_density_matrix_(natoms::Ref{Int32}, atom_rho::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### set_density_matrix_ = ", status[1])

# get compensation charge
@ccall LA.live_atom_get_compensation_charge_(natoms::Ref{Int32}, atom_rho::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### get_compensation_charge_ = ", status[1])

# set potential multipole
@ccall LA.live_atom_set_potential_multipole_(natoms::Ref{Int32}, atom_rho::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### set_potential_multipole_ = ", status[1])

# get zero potential
@ccall LA.live_atom_get_zero_potential_(natoms::Ref{Int32}, atom_rho::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### get_zero_potential_ = ", status[1])
if plot
    for ia = 1:natoms
        zero_potential = atom_rho[ia]
        println("### r, zero_potential(r) for atom #", ia)
        for ir2 = 1:nr2
            println(sqrt((ir2 - 1)/ar2), " ", zero_potential[ir2])
        end # ir2
        println("")
    end # ia
end # plot


function nSHO(numax)
    return Base.div((numax + 1)*(numax + 2)*(numax + 3), 6)
end # nSHO

# get projectors
@ccall LA.live_atom_get_projectors_(natoms::Ref{Int32}, sigma::Ref{Float64}, numax::Ref{Int32}, status::Ref{Int32})::Cvoid
println("### get_projectors_ = ", status[1])
if true
    for ia = 1:natoms
        println("### SHO projectors for atom #", ia," numax= ", numax[ia], " sigma= ", sigma[ia], " Bohr")
    end # ia
end # plot

# get hamiltonian matrix
@ccall LA.live_atom_get_hamiltonian_matrix_(natoms::Ref{Int32}, pointers::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### get_hamiltonian_matrix_ = ", status[1])
if plot
    using Formatting
    for ia = 1:natoms
        nsho = nSHO(numax[ia])
        println("### ",nsho," x ",nsho," Hamiltonian for atom #", ia," in Hartree")
        Hmt = pointers[ia]
        for i = 1:nsho
            print("### i=",i,"  ")
            for j = 1:nsho
                print(Formatting.sprintf1(" %.3f", Hmt[(i - 1)*nsho + j]))
            end # j
            println()
        end # i
    end # ia
end # plot

# get energy contributions
energy = zeros(Float64, natoms)
@ccall LA.live_atom_get_energy_contributions_(natoms::Ref{Int32}, energy::Ref{Float64}, status::Ref{Int32})::Cvoid
println("### get_energy_contributions_ = ", status[1])
if true
    for ia = 1:natoms
        println("### total energy contribution for atom #", ia," is ", energy[ia], " Hartree")
    end # ia
end # true

# direct update (developer access)
ip = zeros(Int32, 1)
dp = zeros(Float64, 1)
fp = zeros(Float32, 1)
@ccall LA.live_atom_update_("direct"::Cstring, natoms::Ref{Int32}, dp::Ref{Float64}, ip::Ref{Int32}, fp::Ref{Float32}, pointers::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### update_ = ", status[1])


@ccall LA.live_atom_finalize_(natoms::Ref{Int32}, status::Ref{Int32})::Cvoid
println("### live_atom_finalize_ = ", status[1])
