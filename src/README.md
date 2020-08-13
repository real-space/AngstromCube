**A4cube**

*This README is only about the code base, please refer to ../README.md for the entire application*
  
AAcube code base is structured as follows:
There is at max one namespace per file.
Each namespace offers the function
```C++
    status_t all_tests(int const echo=1);
```
and should be listed in *main.cxx* to be envoked via the command line using the `-t` option.
The executable shows a simple help when called with `-h`.

Output to *stdout* should always depend on the `echo`-levels so it can be turned off completely.
Warnings printed to stdout, to stderr, recorded in memory and summarized at the programs end 

All log output to *stdout* should start a line with `# ` and, if the information is about 
a specific atom, followed by some atom identifier label.

For development it is sometimes helpful to have plotable output. Here are two options:
  1. Use the helper functions `dump_to_file` offered by *debug_output.hxx*
  2. Print to stdout directly in the format `x y0 y1 y2 ... yn` without `# ` in front
The latter can easily be plotted e.g. by the `xmgrace -nxy` tool.
For clarity, prepend a line with `## x_legend y_legends :` (mind the double hashtag)
When plotting an entire log file `grep \#\# <logfile>` can then help to identify graphs.

When including functionality from another header file (self-written onces or libraries),
the coder is encouraged to keep a list of symbols up to date, listing only those
symbols that re used in the compilation unit at hand.
```C++
    #include "debug_output.hxx" // dump_to_file
    #include "debug_tools.hxx" // ::read_from_file
```
The syntax *::read_from_file* indicates that there is a namespace in front
which matches with the filename of the header file, i.e. *debug_tools* in this example.
This convention is inspired by the *use <module>, only: <symbols>* syntax of Fortran90. 

*Macros and Preprocessor symbols*
The most important preprocessor symbol is
```C
  #ifdef DEVEL
  // code that is useful during development but should not appear in release versions
  #endif
```
The *cppp -U DEVEL* tool is used to remove all DEVEL-code from the code base while
maintaining the other macros. Here, *-U* stands for undefine.
Please try to keep the release code base small and readable!
