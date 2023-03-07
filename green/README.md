*Green function code*

Compilation

    green can be compiled for the CPU without CUDA to check its functionality.
    Otherwise, green is a CUDA code for NVIDIA GPUs.
    Please use either the Makefile in this folder or CMake

Input

    green expects Hmt.xml, an XML-formatted ASCII file containing
    the local and non-local parts of the potential operator as
    exported by write_to_file in include/grid_operators.hxx

Dependencies

    green needs the header-only library rapidxml.
    green can include the header-only library tfQMRgpu.
