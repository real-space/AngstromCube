## Example from https://juliapackages.com/p/lightxml
# import Pkg; Pkg.add("LightXML")
using LightXML

# parse example.xml:
# xdoc is an instance of XMLDocument, which maintains a tree structure
xdoc = parse_file("example.xml")

# get the root element
xroot = root(xdoc)  # an instance of XMLElement
# print its name
# println(name(xroot))  # this should print: bookstore

# traverse all its child nodes and print element names
for c in child_nodes(xroot)  # c is an instance of XMLNode
#     println(nodetype(c))
    if is_elementnode(c)
        e = XMLElement(c)  # this makes an XMLElement instance
#         println(name(e))
    end
end

#=
If the remainder of the script does not use the document or any of its children,
you can call free here to deallocate the memory. The memory will only get
deallocated by calling free or by exiting julia -- i.e., the memory allocated by
libxml2 will not get freed when the julia variable wrapping it goes out of
scope.
=#
free(xdoc)
