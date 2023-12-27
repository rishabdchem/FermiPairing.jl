using LinearAlgebra 
using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export getesp, getpermanent
export getbtp, btpdim

#====================================
Common functions.
====================================#

"""
Compute the scalar corresponding to the elementary 
symmetric polynomial (ESP) of a vector.

For the `SumESP` algorithm, see
DOI: 10.1016/j.amc.2015.08.134.

# Arguments
- `V`: vector of dimension m 
- `n`: degree of ESP 
"""
function getesp(V::AbstractVector{T}, 
                n::Int) where T<:Number 

    # Check
    m = size(V, 1)
    m < 0 && error("m < 0 in getesp")    
    n < 0 && error("n < 0 in getesp")    

    # Trivial
    m < n && return zero(T) 
    n == 0 && return one(T) 

    # Initialize
    E = zeros(T, m, n+1)
    E[:, 1] = ones(T, m) 
    E[1, 2] = V[1]

    # SumESP
    for i = 2:m, j = max(1, i + n - m):min(i, n)
        E[i, j+1] = E[i-1, j+1] 
        E[i, j+1] += V[i] * E[i-1, j]
    end

    return E[m, n+1]
end


"""
Compute the permanent of a square matrix. 
Brute-force computation.

# Argument
- `A`: matrix of n x n dimensions 
"""
function getpermanent(A::AbstractMatrix{T}) where T<:Number 

    # Check
    size(A, 1) == size(A, 2) || error("Wrong matrix dimension in getpermanent")    

    # Initialize
    n = size(A, 1)
    val = zero(T)
    perm = permutations(Vector(1:n)) 

    # Compute
    for mu in perm
        intval = one(T)
        for j = 1:n
            p = mu[j]
            intval *= A[j, p]
        end 
        val += intval 
    end

    return val
end


"""
Solve H C = M C E for the ground state.

H and M must be Hermitian.
M must be positive (semi) definite. 

# Arguments
- `H`: principal matrix 
- `M`: metric matrix 
- `tolmax`: tolerance for declaring near-zero modes
"""
function nocigs(H::AbstractMatrix{T}, 
		M::AbstractMatrix{T}=Matrix{T}(I, size(H, 1), size(H, 1)), 
	        tolmax::Float64=1.0e-12) where T<:BlasFloat

    # Check
    size(H) == size(M) || error("H and M size must match in nocigs")
    size(H, 1) == size(H, 2) || error("H must be a square matrix in nocigs")
    size(H, 1) < 2 && error("H dim < 2 in nocigs")
    ishermitian(H) || error("H is not Hermitian in nocigs")
    ishermitian(M) || error("M is not Hermitian in nocigs")

    # Eigen decomposition of M	  
    evalsm, evecsm = eigen(M) 

    # Count near-zero modes of M
    zeromd, isallnegative, isboth = eigenanalysis(evalsm, tolmax) 

    # Check
    isallnegative == true && error("Check metric sign in nocigs")   
    isboth == true && error("Unphysical metric in nocigs")   
        
    # Choose generalized diag route for H
    if zeromd == 0
         evalsh, evecsh = eigen(H, M) 
    else
         evalsh, evecsh = zeromdsolver(H, evalsm, evecsm, zeromd) 
    end

    # GS results 
    return evalsh[1], evecsh[:, 1], zeromd
end


"""
General H C = M C E.

H and M must be Hermitian.
M must be positive (semi) definite. 

# Arguments
- `H`: principal matrix 
- `M`: metric matrix 
- `nstate`: number of modes we want 
- `tolmax`: tolerance for declaring near-zero modes
"""
function nociall(H::AbstractMatrix{T}, 
		 M::AbstractMatrix{T}=Matrix{T}(I, size(H, 1), size(H, 1)), 
                 nstate::Int=1,
	         tolmax::Float64=1.0e-12) where T<:BlasFloat

    # Check
    size(H) == size(M) || error("H and M size must match in nociall")
    size(H, 1) == size(H, 2) || error("H must be a square matrix in nociall")
    size(H, 1) < 2 && error("H dim < 2 in nociall")
    ishermitian(H) || error("H is not Hermitian in nociall")
    ishermitian(M) || error("M is not Hermitian in nociall")

    # Eigen decomposition of M	  
    evalsm, evecsm = eigen(M) 

    # Count near-zero modes of M
    zeromd, isallnegative, isboth = eigenanalysis(evalsm, tolmax) 
    zeromd + nstate > size(H, 1) && error("Zeromd + nstate > eigen dimension in nociall")

    # Check
    isallnegative == true && error("Check metric sign in nocigs")   
    isboth == true && error("Unphysical metric in nocigs")   
        
    # Choose generalized diag route for H
    if zeromd == 0
         evalsh, evecsh = eigen(H, M) 
    else
         evalsh, evecsh = zeromdsolver(H, evalsm, evecsm, zeromd) 
    end

    # Copy results 
    return evalsh[1:nstate], evecsh[:, 1:nstate], zeromd
end


"""
C <== H C = M C E
given zero mode info on PSD Hermitian metric M.

H and M must be Hermitian.

# Arguments
- `H`: matrix that would be transformed 
- `D`: eigenvalues of metric 
- `U`: eigenvectors of metric 
- `zeromd`: number of zero metric modes 
"""
function zeromdsolver(H::AbstractMatrix{T}, 
		      D::AbstractVector{Float64}, 
		      U::AbstractMatrix{T}, 
		      zeromd::Int) where T<:BlasFloat

    # Check
    zeromd < size(H, 1) - 1 || error("Too many zero modes in zeromdsolver") 
    ishermitian(H) || error("H is not Hermitian in zeromdsolver")

    # Transformed H and coefficient transformer 
    newh, matx = canonortho(H, D, U, zeromd) 

    # Check
    size(newh, 1) == size(H, 1) - zeromd || 
        error("newh dimension error from canonortho") 
    size(matx, 1) == size(H, 1) ||
        error("matx dimension error from canonortho") 
    size(matx, 2) == size(H, 1) - zeromd || 
        error("matx dimension error from canonortho") 
    ishermitian(newh) || error("newh is not Hermitian in zeromdsolver")

    # Eigen decomposition of newh 
    evalsh, evecsh = eigen(newh) 

    # Transform evecsh
    evecsh = matx * evecsh

    return evalsh, evecsh
end


"""
Canonical orthogonalization of matrix given 
its Hermitian PSD metric info. 
See section 3.4.5 of Modern Quantum Chemistry 
by Szabo and Ostlund, Dover (1996).

H and M must be Hermitian.

# Arguments
- `H`: matrix that would be transformed 
- `D`: eigenvalues of metric 
- `U`: eigenvectors of metric 
- `zeromd`: number of zero metric modes 
"""
function canonortho(H::AbstractMatrix{T}, 
		    D::AbstractVector{Float64}, 
		    U::AbstractMatrix{T}, 
		    zeromd::Int) where T <:BlasFloat

    # Initialize
    matx = zeros(T, size(H, 1), size(H, 1) - zeromd)

    # X = U D^{-1/2} considering zero modes
    for q = 1:(size(H, 1) - zeromd), p = 1:size(H, 1)
        matx[p, q] = U[p, q + zeromd] / sqrt(D[q + zeromd])
    end

    # H <= X! H X
    intmat = H * matx
    newh = matx' * intmat

    # Output X and transformed H
    return Hermitian(newh), matx 
end


"""
Eigenvalue analysis (real part only) of a matrix.

Does the matrix have *zero* modes?
Does the matrix have all *negative* modes?
Does the matrix have both *positive* and *negative* modes?

# Arguments
- `D`: eigenvalues of input matrix 
- `tolmax`: tolerance for declaring near-zero modes
"""
function eigenanalysis(D::AbstractVector{T}, 
	               tolmax::Float64=1.0e-12) where T<:BlasFloat 

    zeromd = 0
    indn = 0
    indp = 0
    for p = 1:size(D, 1)
        if abs(D[p]) < tolmax 
            zeromd += 1     
        elseif real(D[p]) < 0.0
            indn += 1 
        elseif real(D[p]) > 0.0
            indp += 1 
        end
    end
   
    isallnegative = true 
    isboth = false 
    if indp + zeromd > 0 
        isallnegative = false 
        if indn > 0 
            isboth = true
        end
    end

    # Return number of zero modes and analysis    
    return zeromd, isallnegative, isboth 
end

"""
Compute the scalar corresponding to the binary 
tree polynomial (BTP) of a matrix.

# Argument
- `A`: matrix of dimensions n times m 
"""
function getbtp(A::AbstractMatrix{T}) where T<:Number 

    # Check
    m = size(A, 2)
    n = size(A, 1)
    m < 0 && error("m < 0 in getbtp")    
    n < 0 && error("n < 0 in getbtp")    

    # Trivial
    m < n && return zero(T) 
    n == 0 && return one(T) 

    # Initialize
    B = zeros(T, m, n+1)
    B[:, 1] = ones(T, m) 
    B[1, 2] = copy(A[1, 1])

    # SumBTP
    for i = 2:m, j = btprowmin(m, n, i):btprowmax(n, i) 
        B[i, j+1] = copy(B[i-1, j+1]) 
        B[i, j+1] += ( A[j, i] * B[i-1, j] )
    end

    return B[m, n+1]
end


"""
Compute all BTP(p, j) of an n x m matrix, where p = 1-m and j = 0-n.  

# Argument
- `A`: matrix of dimensions n times m 
"""
function getbtpmat(A::AbstractMatrix{T}) where T<:Number 

    # Check
    m = size(A, 2)
    n = size(A, 1)

    # Initialize
    B = zeros(T, m, n+1)
    B[:, 1] = ones(T, m) 
    B[1, 2] = copy(A[1, 1])

    # SumBTP
    for i = 2:m, j = 1:n 
        i < j && continue
        B[i, j+1] = copy(B[i-1, j+1]) 
        B[i, j+1] += ( A[j, i] * B[i-1, j] )
    end

    return B
end


"""
Unpack BTP vector followed by any permutation.

# Arguments
- `V`: BTP vector 
- `m`: order of BTP 
- `n`: degree of BTP 
- `indper`: permutation of rows
"""
function btpvectomat(V::AbstractVector{T},
                     m::Int,
                     n::Int, 
                     indper::Vector{Int}=Vector(1:n)) where T<:Number

    # Check
    d = size(V, 1)
    m == 0 && error("Wrong m in btpvectomat")
    n == 0 && error("Wrong n in btpvectomat")
    d == btpdim(m, n) || error("Wrong vector dimension in btpvectomat")

    # Initialize
    G = zeros(T, n, m)
    X = zeros(T, n, m)
    G[1, 1] = V[1]

    # Unpack
    d = 1
    for i = 2:m, j = btprowmin(m, n, i):btprowmax(n, i) 
        d += 1     
        G[j, i] = V[d]
    end
    indper == Vector(1:n) && return G

    # Permute 
    for j = 1:n
        X[indper[j], :] .= copy(G[j, :])
    end
   
    return X 
end


"""
Pack BTP matrix after "de-permutation".

# Argument
- `G`: BTP matrix of dimensions n times m 
- `indper`: permutation of rows
"""
function btpmattovec(G::AbstractMatrix{T},
                     indper::Vector{Int}=Vector(1:size(G, 1))) where T<:Number

    # Check 
    m = size(G, 2)
    n = size(G, 1)
    m == 0 && error("Wrong m in btpmattovec")
    n == 0 && error("Wrong n in btpmattovec")
    d = btpdim(m, n)

    # Initialize
    V = zeros(T, d)
    X = zeros(T, n, m)

    # Permute 
    if indper != Vector(1:n) 
        for j = 1:n
            X[j, :] .= copy(G[indper[j], :])
        end
        G .= copy(X)
    end 

    # Pack
    V[1] = G[1, 1] 
    d = 1
    for i = 2:m, j = btprowmin(m, n, i):btprowmax(n, i) 
        d += 1     
        V[d] = G[j, i] 
    end
   
    return V 
end


"""
Number of unique BTP variables.

# Arguments
- `m`: order of BTP 
- `n`: degree of BTP 
"""
function btpdim(m::Int,
                n::Int)

    d = 1
    for i = 2:m, j = btprowmin(m, n, i):btprowmax(n, i) 
        d += 1     
    end
   
    return d
end


function btprowmin(m::Int,
                   n::Int,
                   p::Int)

    return max(1, p + n - m)
end


function btprowmax(n::Int,
                   p::Int)

    return min(p, n)
end


"""
C2 rotation of a rectangular matrix perpendicular to its plane.

# Argument
- `A`: rectangular matrix input 
"""
function matc2rot(A::AbstractMatrix{T}) where T<:Number 

    B = zeros(T, size(A, 1), size(A, 2))
    for j = 1:size(A, 1), p = 1:size(A, 2) 
        jmod = size(A, 1) - j +1
        pmod = size(A, 2) - p +1
        B[jmod, pmod] = copy(A[j, p])
    end

    return B
end


