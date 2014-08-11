# Interface to the Watson Sparse Matrix Package from IBM

The
[Watson Sparse Matrix Package](http://research.ibm.com/projects/wsmp)
is proprietary software. Use of this [Julia](http://julialang.org/)
interface requires that the library `libwsmp64.a` be downloaded and
installed in the directory
```julia
Pkg.dir("WSMP","deps","usr","lib")
```
and that a valid `wsmp.lic` license file be available in the same
directory.

The library is shipped as a static library.  It must be converted to a
dynamic library before use. A rudimentary `build.jl` file is available
to accomplish this.

Note that `libwsmp.so` __must__ be linked against a BLAS compiled with
the environment variable `USE_BLAS` set to `0`. An easy way to
accomplish this is to add the line
```
USE_BLAS64=0
```
in the `Make.user` file in the Julia home directory and run
```
make cleanall
make
make testall
```

Check that BLAS are compiled to use 32-bit integers with
```
julia> Base.LinAlg.BlasInt
Int32
```

## WSMP types

The `Wssmp` type represents a symmetric sparse matrix. A constructor
that takes a symmetric `SparseMatrixCSC` is available.

By default the Wssmp type uses a CSC (compressed sparse column)
representation of the lower triangle of the original matrix.
The MSC (modified compressed sparse column) representation, in which
the diagonal is stored separately from the strict lower triangle in
CSC format, is also available.

