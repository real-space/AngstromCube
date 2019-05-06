How can a method based on Cartesian real-space grids treat arbitrary unit cells?

input Bravais matrix
    a11 a12 a13
    a21 a22 a23
    a31 a32 a33

Assert that a11 is non-zero
Always transform into an upper triangular matrix by row ops
    a11 a12 a13
     0  b22 b23
     0   0  c33

    b22= a22 - a21*a12/a11
    b23= a23 - a21*a13/a11

    b32= a32 - a31*a12/a11
    b33= a33 - a31*a13/a11

Assert that b22 is non-zero

    c33= b33 - b32*b23/b22

Then, we have a geometry which can be represented by a rectangular cell (a11 x b22 x c33)
with shift-boundary-conditions according to a12, a13 and b23

