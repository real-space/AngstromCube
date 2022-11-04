######################################################
### interface the LiveAtom library
######################################################

const LA = "libliveatom.so"

plot = false

status = 0
@ccall LA.live_atom_set_env_("called.from"::Cstring, "Julia"::Cstring, status::Ref{Int32})::Cvoid
# ccall(("live_atom_set_env_", LA), Cvoid, (Cstring, Cstring, Ref{Int32}), "called.from", "Julia", status)
println("### live_atom_init_env = ", status)

status = 0
@ccall LA.live_atom_init_env_("control.sh"::Cstring, status::Ref{Int32})::Cvoid
# ccall(("live_atom_init_env_", LA), Cvoid, (Cstring, Ref{Int32}), "control.sh", status)
println("### live_atom_init_env_ = ", status)

# the total number of atoms
const natoms = 1

# initialize the atoms using pot/Zeff.00Z files
Z_core = ones(Float64, natoms)
atom_id = zeros(Int32, natoms)
numax = zeros(Int32, natoms)
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
status = 0
for ia = 1:natoms
    Z_core[ia] = 10 + ia
    atom_id[ia] = ia - 1
end # ia
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
println("### live_atom_initialize_ = ", status)


# get the smooth core density
const ar2 = 16.
const nr2 = 2^12
pointers = Vector{Vector{Float64}}(undef, natoms)
for ia = 1:natoms
    pointers[ia] = zeros(Float64, nr2)
end # ia
status = 0
@ccall LA.live_atom_get_core_density_(natoms::Ref{Int32}, pointers::Ref{Ptr{Float64}}, status::Ref{Int32})::Cvoid
println("### live_atom_get_core_density_ = ", status)
if plot
    for ia = 1:natoms
        core_density = pointers[ia]
        println("### r, core_density(r) for atom #", ia)
        for ir2 = 1:nr2
            println(sqrt((ir2 - 1)/ar2), " ", core_density[ir2])
        end # ir2
        println("")
    end # ia
    println("")
end # plot



status = 0
@ccall LA.live_atom_finalize_(natoms::Ref{Int32}, status::Ref{Int32})::Cvoid
println("### live_atom_finalize_ = ", status)
