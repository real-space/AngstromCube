using LightXML

### ToDo: move to sho_tools
nSHO(numax::Integer) = Base.div((numax + 1)*(numax + 2)*(numax + 3), 6)

convert(str::AbstractString) = Base.parse(Float64, str)
average(array) = Base.sum(array)/Base.length(array)

echo = false

# parse example.xml:
# xdoc is an instance of XMLDocument, which maintains a tree structure
filename = "Hmt.xml"
xdoc = LightXML.parse_file(filename)

# get the root element
grid_Hamiltonian = LightXML.root(xdoc)  # an instance of XMLElement
println("# ",LightXML.name(grid_Hamiltonian)," version=",LightXML.attribute(grid_Hamiltonian, "version"))  # this should print: grid_Hamiltonian

sho_atoms = LightXML.find_element(grid_Hamiltonian, "sho_atoms")
atoms = LightXML.get_elements_by_tagname(sho_atoms, "atom")
natoms = Base.parse(Int64, LightXML.attribute(sho_atoms, "number"))
# println("# length(atoms)=",Base.length(atoms)," natoms=",natoms)
@assert Base.length(atoms) == natoms "XML file inconsistent"
apos = Base.zeros(Float64, 8, natoms)
for atom in atoms
    gid = Base.parse(Int64, LightXML.attribute(atom, "gid"))
    println("# Atom gid=",gid)
    @assert gid < natoms
    position = LightXML.find_element(atom, "position")
    apos[1:3,gid + 1] = convert.(LightXML.attribute.([position, position, position], ["x", "y", "z"]))
    projectors = find_element(atom, "projectors")
    @assert "sho" == LightXML.attribute(projectors, "type")
    numax = Base.parse(Int8, LightXML.attribute(projectors, "numax"))
    apos[5,gid + 1] = gid
    apos[6,gid + 1] = numax
    apos[7,gid + 1] = convert(LightXML.attribute(projectors, "sigma"))
    nsho = nSHO(numax)
    println("# numax= ",numax," expects ",nsho," coefficients")
    hamiltonian = LightXML.find_element(atom, "hamiltonian")
    hamiltonian_values = Base.split(LightXML.content(hamiltonian))
    @assert Base.length(hamiltonian_values) == nsho^2
    hamiltonian_values = Base.reshape(convert.(hamiltonian_values), (nsho, nsho))
    overlap = LightXML.find_element(atom, "overlap")
    overlap_values = Base.split(LightXML.content(overlap))
    @assert Base.length(overlap_values) == nsho^2
    overlap_values = Base.reshape(convert.(overlap_values), (nsho, nsho))
    if echo
        println("# Hamiltonian values = ",hamiltonian_values)
        println("# Overlap values     = ",overlap_values)
    end
end # atom
println("# xyzZinso = ",apos)

spacing = LightXML.find_element(grid_Hamiltonian, "spacing")
grid_spacing = convert.(LightXML.attribute.([spacing, spacing, spacing], ["x", "y", "z"]))
println("# grid spacing = ",grid_spacing," Bohr")

potential = LightXML.find_element(grid_Hamiltonian, "potential")
ng = Base.parse.(Int64, LightXML.attribute.([potential, potential, potential], ["nx", "ny", "nz"]))
println("# Potential shape is ",ng," grid points")
potential_values = Base.split(LightXML.content(potential))
@assert Base.length(potential_values) == ng[1]*ng[2]*ng[3]
# println(potential_values) ## potentially long output
potential_values = Base.reshape(convert.(potential_values), (ng[1], ng[2], ng[3])) # apply convert elementwise
# println(potential_values) ## potentially long output
println("# Average potential value = ",average(potential_values)," Hartree")

free(xdoc)
