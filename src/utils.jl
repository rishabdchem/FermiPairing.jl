using LinearAlgebra 

import LinearAlgebra.BLAS: BlasFloat

export getesp

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
    size(H, 1) == size(H, 2) || error("H must be a square matrix")
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

    # Return GS eval, evec, and number of zero modes
    return evalsh[1], evecsh[:, 1], zeromd
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



