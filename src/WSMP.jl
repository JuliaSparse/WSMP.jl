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

include("wssmp.jl")

end # module
