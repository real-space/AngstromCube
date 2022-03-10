*Green function code*

Compilation

    green can be compiled for the CPU without CUDA to check the functionality.
    To do this, create a soft link
        ln -s green.cu green.cxx
    
    Otherwise, green is a CUDA code for NVIDIA GPUs.
    
Input

    green expects Hmt.xml, an XML-formatted ASCII file containing
    the local and non-local parts of the potential operator.

Dependencies

    currently green needs the header-only library rapidxml.
    planned extension foresees the inclusion of tfQMRgpu.
