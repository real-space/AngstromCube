import LightXML: name, attribute, find_element, get_elements_by_tagname, parse_file, root, content, free

### ToDo: move to sho_tools
nSHO(numax::Integer) = div((numax + 1)*(numax + 2)*(numax + 3), 6) # number of SphericalHarmonicOscillator functions up to numax

convert(str::AbstractString) = parse(Float64, str)
average(array) = sum(array)/length(array)

echo = false

# parse example.xml:
# xdoc is an instance of XMLDocument, which maintains a tree structure
filename = "Hmt.xml"
xdoc = parse_file(filename)

# get the root element
grid_Hamiltonian = root(xdoc)  # an instance of XMLElement
println("# ",name(grid_Hamiltonian)," version=",attribute(grid_Hamiltonian, "version"))  # this should print: grid_Hamiltonian

sho_atoms = find_element(grid_Hamiltonian, "sho_atoms")
natoms = parse(Int64, attribute(sho_atoms, "number"))
atoms = get_elements_by_tagname(sho_atoms, "atom")
# println("# length(atoms)=",length(atoms)," natoms=",natoms)
@assert length(atoms) == natoms
apos = zeros(Float64, 8, natoms)
for atom in atoms
    gid = parse(Int64, attribute(atom, "gid"))
    println("# Atom gid=",gid)
    @assert gid < natoms
    position = find_element(atom, "position")
  # apos[1:3,gid + 1] = convert.(attribute.(Ref(position), ["x", "y", "z"]))
    for (i,xyz) in zip(1:3, ["x", "y", "z"])
        apos[i,gid + 1] = convert(attribute(position, xyz))
    end
    projectors = find_element(atom, "projectors")
    @assert "sho" == attribute(projectors, "type")
    numax = parse(Int32, attribute(projectors, "numax"))
    apos[5,gid + 1] = gid
    apos[6,gid + 1] = numax
    apos[7,gid + 1] = convert(attribute(projectors, "sigma"))
    nsho = nSHO(numax)
    println("# numax= ",numax," expects ",nsho," coefficients")
    hamiltonian = find_element(atom, "hamiltonian")
    hamiltonian_values = split(content(hamiltonian))
    @assert length(hamiltonian_values) == nsho^2
    hamiltonian_values = reshape(convert.(hamiltonian_values), (nsho, nsho))
    overlap = find_element(atom, "overlap")
    overlap_values = split(content(overlap))
    @assert length(overlap_values) == nsho^2
    overlap_values = reshape(convert.(overlap_values), (nsho, nsho))
    if echo
        println("# Hamiltonian values = ",hamiltonian_values)
        println("# Overlap values     = ",overlap_values)
    end # echo
end # atom
println("# xyzZinso = ",apos)

spacing = find_element(grid_Hamiltonian, "spacing")
# grid_spacing = convert.(attribute.(Ref(spacing), ["x", "y", "z"]))
grid_spacing = [1.0, 1.0, 1.0]
for (i,xyz) in zip(1:3, ["x", "y", "z"])
    grid_spacing[i] = convert(attribute(spacing, xyz))
end
println("# grid spacing = ",grid_spacing," Bohr")

potential = find_element(grid_Hamiltonian, "potential")
# ng = parse.(Int64, attribute.(Ref(potential), ["nx", "ny", "nz"]))
ng = [1, 1, 1]
for (i,xyz) in zip(1:3, ["nx", "ny", "nz"])
    ng[i] = parse(Int64, attribute(potential, xyz))
end
println("# Potential shape is ",ng," grid points")
potential_values = split(content(potential))
@assert length(potential_values) == ng[1]*ng[2]*ng[3]
# println(potential_values) ## potentially long output
potential_values = reshape(convert.(potential_values), (ng[1], ng[2], ng[3]))
# println(potential_values) ## potentially long output
println("# Average potential value = ",average(potential_values)," Hartree")

free(xdoc)
