***Multi-dimensional Arrays***

The largest feature-envy for Fortran90, however, are the intrinsic support of
multi-dimensional arrays with support for run-time checking of index-out-of-bounds.

Therefore, `data_view.hxx` introduces the containers
```C++
    view2D<T> , view3D<T>, /* and */ view4D<T>
```
These can be initialized similar to `std::vector<T>` using allocating constructors
described in the following.

*Allocating constructors*
```C++
    view2D<T> a2(n1, n0);         // allocates       n1*n0*sizeof(T) Bytes
    view3D<T> a3(n2, n1, n0);     // allocates    n2*n1*n0*sizeof(T) Bytes
    // and
    view4D<T> a4(n3, n2, n1, n0); // allocates n3*n2*n1*n0*sizeof(T) Bytes
```
and initializer values can be passed, e.g.
```C++
    view2D<float> a2(n1, n0, 1.0f); // allocates n1*n0 floats, set all to value 1
```
We can use the familiar `[]`-operator onto a `view2D`.
It does not return an `std::vector` or a `view1D` for reasons of performance
but a plain `T*` pointer where it is the programmer`s duty to check that
this pointer is not dereferenced beyond entry `n0`, i.e. index checking is disabled here.

The `[]`-operator onto a `view3D` returns a `view2D` and, similarly,
for `view4D` a `view3D` is retained. Mind that the returned objects 
are only data views (hence the name) and do not allocate their own memory.
This allows to use, if necessary, also the `[]`-indexing syntax in a cascading fashion
```C++
    view4D<T> a4(n3, n2, n1, n0); // allocates
    T b = a4[i3][i2][i1][i0];
```
Although deep copies are avoided, it comes at a runtime overhead of creating 
several temporary objects. It is encouraged to use the Fortran-like `()`-indexing syntax
```C++
    b = a4(i3,i2,i1,i0);
```
to minimize the overhead.
Furthermore, the `()`-operator allows to pass 1D sub-views of e.g. 3D arrays
```C++
    view3D<double> a3(n2, n1, n0); // allocates memory
    double* a1 = a3(i2, i1);
```
or 2D sub-views of 4D arrays as
```C++
    view4D<double> a4(n3, n2, n1, n0); // allocates memory
    view2D<double> a2 = a4(i3, i2);
```
The `()`-operator with a single index is equivalent to the `[]`-operator mentioned above.

*Alternative constructors*
As described for the `[]`-operator for `view3D` and `view4D`,
we retrieve only views to the data for which the index computation is known.
An example
```C++
    view4D<double> a3(n2, n1, n0, 5.0); // allocates memory and sets all to 5.0
    view2D<double> a2 = a3[i2];
    a2(0,0) = 3.0;
    return a3(i2, 0, 0);
```
will return `3.0` although it looks like we only modified a matrix element of `a2`.
Internally, however, `a2` points to overlapping memory regions with `a3`.
In this spirit, there is a second set of constructor methods:

*Wrapping constructors*
```C++
    T* p = new T(n3*n2*n1*n0); // get memory elsewhere
    view2D<T> a2(p, n0);         // wraps pointer p, assumed shape a2(?, n0)
    view3D<T> a3(p, n1, n0);     // wraps pointer p, assumed shape a2(?, n1, n0)
    // and
    view4D<T> a4(p, n2, n1, n0); // wraps pointer p, assumed shape a4(?, n2, n1, n0)
```
As indicated by `?` index-out-of-bounds cannot be checked for the highest dimension
for views that were constructed as a wrap around a plain pointer `p`.
Also, the `view`s are not memory owners, i.e. they are not responsible for
freeing the memory when destructed.

When passing `view`s to external interfaces, we can access the memory directly by
```C++
    view2D<double> a2(N, M); // allocates N*M doubles
```

