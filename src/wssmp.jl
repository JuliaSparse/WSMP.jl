type Wssmp                              # WSMP symmetric type
    ia::Vector{Cint}
    ja::Vector{Cint}
    avals::Vector{Cdouble}
    diag::Vector{Cdouble}
    perm::Vector{Cint}
    invp::Vector{Cint}
    mrp::Vector{Cint}
    iparm::Vector{Cint}
    dparm::Vector{Cdouble}
end

function Wssmp(ia::Vector{Cint},ja::Vector{Cint},avals::Vector{Cdouble},
                 dd::Vector{Cdouble}=Cdouble[])
    1 == ia[1] || error("Column pointers must be 1-based.  That is, ia[1] == 1")
    (n = length(ia) - 1) > 0 || error("0 × 0 sparse matrices not allowed")
    (ld = length(dd)) == 0 || ld == n ||
        error("length(diag) = $length(diag), should be 0 or $n")
    length(ja) == length(avals) == ia[end] - 1 || throw(DimensionMismatch(""))
    for j in ja
        0 < j ≤ n || error("elements of ja must be in [1,$n]")
    end
    sorted = true
    for i in 1:n         # check for sorted columns and lower triangle
        ia[i] ≤ ia[i+1] || error("elements in ia must be non-decreasing")
        for k in ia[i]:(ia[i+1]-1)
            (ld == 0 ? i ≤ ja[k] : i < ja[k]) || error("matrix is not lower triangular")  
            k == ia[i] || ja[k-1] ≤ ja[k] || (sorted = false)
        end
    end
    if !sorted                          # sort columns if necessary
        info = zeros(Cint,1)
        ccall((:ws_sortindices_d__,libwsmp),Void,
              (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cint}),
              &n,&n,ia,ja,avals,info)
        info[1] == 0 || error("error code $(info[1]) from ws_sortindices_d__")
    end
    ans = Wssmp(ia,ja,avals,dd,Array(Cint,(n,)),Array(Cint,(n,)),Cint[],
                  zeros(Cint,64),zeros(64))
    wssmp(ans,0)
    ld > 0 && (ans.iparm[4] = 1)
    ans
end

function Wssmp(S::SparseMatrixCSC{Cdouble,Cint}, MSC::Bool=false)
    issym(S) || error("Matrix S is not symmetric")
    A = MSC ? tril(S,-1) : tril(S)
    Wssmp(A.colptr,A.rowval,A.nzval,MSC ? diag(S) : Cdouble[])
end
    
function wssmp(a::Wssmp, task::Integer)
    0 ≤ task ≤ 3 || error("task = $task should be in the range [0,3]")
    a.iparm[2] ≤ task || error("wssmp called with task = $task but a.iparm[2] = $(a.iparm[2])")
    a.iparm[3] = task
    task ≡ 0 && (a.iparm[1] = 0)        # basic initialization
    task ≡ 3 && a.iparm[32] ≠ 0 && length(a.diag) ≡ 0 && (a.diag = zeros(length(a.perm)))
    ccall((:wssmp_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},
           Ptr{Cint},Ptr{Cint},Ptr{Void},Ptr{Cint},Ptr{Cint},Ptr{Cint}, Ptr{Cdouble}),
          &(length(a.ia)-1),a.ia,a.ja,a.avals,a.diag,a.perm,a.invp,C_NULL,&1,&0,C_NULL,&0,
          a.mrp,a.iparm,a.dparm)
    if a.iparm[64] ≠ 0
        error("error code $(a.iparm[64]) from ccall to wssmp_ with a.iparm[3] == $(a.iparm[3])")
    end
    a
end

function wssmp(a::Wssmp, b::StridedVecOrMat{Cdouble}, task::Integer=4)
    4 ≤ task ≤ 5 || error("task = $task should be 4 or 5")
    a.iparm[3] = task
    ccall((:wssmp_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},
           Ptr{Cint},Ptr{Cint},Ptr{Void},Ptr{Cint},Ptr{Cint},Ptr{Cint}, Ptr{Cdouble}),
          &(length(a.ia)-1),a.ia,a.ja,a.avals,a.diag,a.perm,a.invp,
          b,&stride(b,2),&size(b,2),C_NULL,&0,a.mrp,a.iparm,a.dparm)
    if a.iparm[64] ≠ 0
        error("error code $(a.iparm[64]) from ccall to wssmp_ with a.iparm[3] == $(a.iparm[3])")
    end
    b
end

Base.A_ldiv_B!(A::Wssmp, b::StridedVecOrMat{Cdouble}) = (A.iparm[2] = min(A.iparm[2],4);wssmp(A,b,4))
(\)(A::Wssmp, b::StridedVecOrMat{Cdouble}) = (A.iparm[2] = min(A.iparm[2],4);wssmp(A,copy(b),4))

Base.size(A::Wssmp) = (n = length(A.ia) - 1; (n,n))

function Base.size(A::Wssmp,i::Integer)
    i ≤ 0 && throw(BoundsError())
    n = length(A.ia) - 1
    1 ≤ i ≤ 2 && return n
    return 1
end

function *(A::Wssmp,x::Vector{Cdouble})
    A.iparm[4] == 0 || error("Wssmp object A must be in CSC format, not MSC")
    wssmatvec(A.ia,A.ja,A.avals,x,similar(x))
end

## for (nm,sym,op5) in ((:wmmrb,:wmmrb_,3),
##                      (:wkktord,:wkktord_,0))
##     @eval begin
##         function $nm(A::SparseMatrixCSC{Cdouble,Cint})
##             issym(A) || error("Matrix A must be symmetric")
##             adj = Base.SparseMatrix.fkeep!(copy(A), (i,j,x,other)->(i!=j), None)
##             n = size(adj,1)
##             options = zeros(Cint,5)
##             options[1] = 3
##             options[3] = 1
##             options[5] = $op5
##             perm = Array(Cint,n)
##             invp = Array(Cint,n)
##             ccall(($(string(nm,"_")),libwsmp),Void,(Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},
##                                                  Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
##                   &n,adj.colptr,adj.rowval,options,&1,perm,invp,C_NULL,&0)
##             perm,invp
##         end
##     end
## end

function wsmmrb(A::SparseMatrixCSC{Cdouble,Cint},options::Vector{Cint}=Int32[3,0,1,0,3])
    issym(A) || error("Matrix A must be symmetric")
    5 ≤ length(options) || error("options vector must be of length 5 or more")
    adj = Base.SparseMatrix.fkeep!(copy(A), (i,j,x,other)->(i ≠ j), None)
    n = size(adj,1)
    perm = Array(Cint,n)
    invp = Array(Cint,n)
    ccall((:wmmrb_,libwsmp),Void,(Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},
                                  Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
          &n,adj.colptr,adj.rowval,options,&1,perm,invp,C_NULL,&0)
    perm,invp
end

function wskktord(A::SparseMatrixCSC{Cdouble,Cint},options::Vector{Cint}=Int32[3,0,1,0,0])
    issym(A) || error("Matrix A must be symmetric")
    5 ≤ length(options) || error("options vector must be of length 5 or more")
    adj = Base.SparseMatrix.fkeep!(copy(A), (i,j,x,other)->(i ≠ j), None)
    n = size(adj,1)
    perm = Array(Cint,n)
    invp = Array(Cint,n)
    ccall((:wkktord_,libwsmp),Void,(Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},
                                  Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
          &n,adj.colptr,adj.rowval,options,&1,perm,invp,C_NULL,&0)
    perm,invp
end
