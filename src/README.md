**AAcube**
    
AAcube code base is structured as follows:
There is at max one namespace per file.
Each namespace offers the function
```C++
    status_t all_tests(int const echo=1);
```
and should be listed in *main.cxx* to be envoked via the command line using the `-t` option.
The executable shows a simple help when called with `-h`.

Output to *stdout* should always depend on the `echo`-levels so it can be turned off completely.
Warnings are recorded and summarized at the programs end (and when dying from fatal errors, ToDo)
All log output to *stdout* should start a line with `# ` and, if the information is about 
a specific atom, followed by some atom identifier label.

For development it is sometimes helpful to have plotable output. Here are two options:
  1. Use the helper functions `dump_to_file` offered by *debug_output.hxx*
  2. Print to stdout directly in the format `x y0 y1 y2 ... yn` with `# `
The latter can easily be plotted e.g. by the `xmgrace -nxy` tool.
For clarity, prepend a line with `## explanations about x and y:`
When plotting an entire log file `grep \#\# log` can then help to identify graphs. 

