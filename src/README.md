**AngstromCube**

*This README is only about the code base, please refer to ../README.md for the entire application*
  
The AngstromCube code base is structured as follows:
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
symbols that appear in the current source file
```C++
    #include "debug_output.hxx" // dump_to_file
    #include "debug_tools.hxx" // ::read_from_file
```
The syntax *::read_from_file* indicates that there is a namespace in front
which should match with the filename of the header file, i.e. *debug_tools* in this example.
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

*You can write Fortran in any language*
The saying above moques about the *bad* programming style many Fortran programs
are written in. However, there are some nice things in Fortran90 that come handy
so they were re-introduced into the C++ of this project:
```Fortran90
    vec(:) = b(:)
```
However, to keep it as basic as possible, it can be done by
```C++
    #include "inline_math.hxx" // set
    set(vec, n, b);
```
but we have to know the length of the arrays.
Certainly, `set` is templated, so where `x` and `b` are pointers, not necessaily of the same type.
`set` also supports scaling and assignment of a scalar value
```C++
    set(vec, n, vector, scalar_factor);
    set(vec, n, scalar_value);
```
More vector functionality can be found in `inline_math.hxx` as e.g.
```C++
    scale(vec, n, scalar_factor);
    scale(vec, n, vector[, scalar_factor]);
    product(vec, n, vector, vector_factors [, vector_factors][, scalar_factor]);
    add_product(vec, n, vector, scalar_factor);
    add_product(vec, n, vector, vector_factors[, scalar_factor]);
```

The largest feature-envy for Fortran90, however, are the intrinsic support of
multi-dimensional arrays including run-time checking of index-out-of-bounds.
Please see `data_view.md` for more details.

*Style Comment*
In each code file (*.hxx, *.cxx) we try to make dependencies easily identifyable.
In Fortran, this would be a "use <modulename>, only: <function1>, <function2>" statement.
In C++ we treat it as a comment after the included header file
    #include <cstdio> // std::printf, ::snprintf
    #include "tools.hxx" // ::foo, ::bar
We only mention the namespace name when it deviates from the header's file name,
e.g. std:: should be mentioned once here.
