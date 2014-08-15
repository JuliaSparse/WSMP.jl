## Julia wrappers of miscellaneous functions described in section 10

function wssortindicesi(m::Integer,n::Integer,ia::Vector{Cint},ja::Vector{Cint})
    info = Array(Cint,1)
    ccall((:ws_sortindices_i_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
          &m,&n,ia,ja,info)
    info[1] == 0 || error("Error return $(info[1]) from ws_sortindices_i_")
    nothing
end

function wssortindicesd(m::Integer,n::Integer,ia::Vector{Cint},ja::Vector{Cint},
                          avals::Vector{Cdouble})
    info = Array(Cint,1)
    ccall((:ws_sortindices_d_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cint}),
          &m,&n,ia,ja,avals,info)
    info[1] == 0 || error("Error return $(info[1]) from ws_sortindices_i_")
    nothing
end

wssetmaxthreads(n::Integer) = ccall((:wsetmaxthrds_,libwsmp),Void,(Ptr{Cint},),&n)

wssetmaxstack(fstk::Cdouble) = ccall((:wsetmaxstack_,libwsmp),Void,(Ptr{Cdouble},),&fstk)

wssetlf(dlf::Cdouble) = ccall((:wsetlf_,libwsmp),Void,(Ptr{Cdouble},),&dlf)

function wsversion()
    V = Array(Cint,1)
    R = Array(Cint,1)
    M = Array(Cint,1)
    ccall((:wsmp_version_,libwsmp),Void,(Ptr{Cint},Ptr{Cint},Ptr{Cint}),V,R,M)
    (V[1],R[1],M[1])
end

wsclear() = ccall((:wsmp_clear_,libwsmp),Void,())

wsffree() = ccall((:wsffree_,libwsmp),Void,())

wsafree() = ccall((:wsafree_,libwsmp),Void,())

wssfree() = ccall((:wssfree_,libwsmp),Void,())

function wssmatvec(ia::Vector{Cint},ja::Vector{Cint},
                   avals::Vector{Cdouble},x::Vector{Cdouble},b::Vector{Cdouble})
    n = length(ia) - 1
    length(x) == length(b) == n || throw(DimensionMismatch(""))
    length(ja) == length(avals) == (ia[end] - 1) || throw(DimensionMismatch(""))    
    ierr = Array(Cint,1)
    ccall((:wssmatvec_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint}),
          &n,ia,ja,avals,x,b,ierr)
    ierr[1] == 0 || error("Error return $(ierr[1]) from wssmatvec_")
    b
end

function wssmatvec(ia::Vector{Cint},ja::Vector{Cint},avals::Vector{Cdouble},x::Vector{Cdouble})
    wssmatvec(ia,ja,avals,x,similar(x))
end

function wsstoremat(id::Integer)
    0 ≤ id < 64 || error("id = $id should be in 0:63")
    info = Array(Cint,1)
    ccall((:wstoremat_,libwsmp),Void,(Ptr{Cint},Ptr{Cint}),&id,info)
    info[1] == 0 && return id
    info[1] == -210 && error("Slot $id is occupied")
    error("Error code $(info[1]) from wstoremat_")
end

function wsrecallmat(id::Integer)
    0 ≤ id < 64 || error("id = $id should be in 0:63")
    info = Array(Cint,1)
    ccall((:wrecallmat_,libwsmp),Void,(Ptr{Cint},Ptr{Cint}),&id,info)
    info[1] == 0 && return id
    info[1] == -210 && error("Slot $id is empty")
    error("Error code $(info[1]) from wstoremat_")
end
