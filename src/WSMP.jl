module WSMP

## Functions to interface with the Watson Sparse Matrix Package, WSMP

## Note that WSMP is not free software.  An evaluation version is available, see
##  http://researcher.watson.ibm.com/projects/wsmp
## Production use requires the purchase of a license.

if isfile(joinpath(Pkg.dir("WSMP"),"deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("WSMP not properly installed. Please run Pkg.build(\"WSMP\")")
end

ENV["MALLOC_TRIM_THRESHOLD_"] = -1
ENV["MALLOC_MMAP_MAX_"] = 0
ENV["WSMPLICPATH"] = Pkg.dir("WSMP","deps","usr","lib")
blas_set_num_threads(1)

export Wssmp,
    wsafree,                   # release the internal copy of the system matrix
    wsclear,                   # release all internal storage
    wsffree,                   # release the internal numeric factor storage
    wsrecallmat,               # recall saved matrix k as active storage
    wssetlf,                   # set the load factor (0. < lf < 1.)
    wssetmaxstack,             # set the max stack size (relative to 1 Mb stack)
    wssetmaxthreads,           # set maximum number of threads to spawn
    wssfree,                   # release the internal symbolic factor storage
    wssmp,                     # call wssmp_ on matrix A for stage k
    wssortindicesd,            # sort row indices carrying along the non-zero values
    wssortindicesi,            # sort row indices
    wsstoremat,                # store the current active matrix in slot k (0 ≤ k ≤ 63)
    wsversion                  # report the version of WSMP in use (returns 3 Cint's)

include("misc.jl")
include("wssmp.jl")

end # module
