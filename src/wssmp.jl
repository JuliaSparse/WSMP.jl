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
    task == 0 && (a.iparm[1] = 0)
    ccall((:wssmp_,libwsmp),Void,
          (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},
           Ptr{Cint},Ptr{Cint},Ptr{Void},Ptr{Cint},Ptr{Cint},Ptr{Cint}, Ptr{Cdouble}),
          &(length(a.ia)-1),a.ia,a.ja,a.avals,a.diag,a.perm,a.invp,C_NULL,&1,&0,C_NULL,&0,
          a.mrp,a.iparm,a.dparm)
    if a.iparm[64] != 0
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
    if a.iparm[64] != 0
        error("error code $(a.iparm[64]) from ccall to wssmp_ with a.iparm[3] == $(a.iparm[3])")
    end
    b
end
