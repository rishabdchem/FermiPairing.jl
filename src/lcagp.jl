using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export pairinglcagp_gen, pairinglcagp_pivot 

#===========================================
LC-AGP with fixed AGPs.
===========================================#

"""
LC-AGP ground state.

# Arguments
- `etamat`: rows of input AGPs 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairinglcagp_gen(etamat::AbstractMatrix{T},
                          n::Int, 
                          hamp::Float64,
                          hamtype::Symbol=:RBCS,
                          bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    m = size(etamat, 2)
    d = size(etamat, 1)
    M = zeros(T, d, d)
    H = zeros(T, d, d)
    eta = copy(etamat)

    # Normalize AGPs
    for i = 1:d
        eta[i, :] = normalizeagp!(eta[i, :], n)
    end

    # Matrices 
    for i = 1:d
	M[i, i] = agpoverlap(eta[i, :], eta[i, :], n)
	H[i, i] = agphamoverlap(eta[i, :], eta[i, :], n, hamp, hamtype, bound)
    end	    
    for i = 1:d-1, j = i+1:d
	M[i, j] = agpoverlap(eta[i, :], eta[j, :], n)
        M[j, i] = M[i, j] 
        H[i, j] = agphamoverlap(eta[i, :], eta[j, :], n, hamp, hamtype, bound)
        H[j, i] = H[i, j] 
    end	    

    # Solve only for ground state 
    @time eval, evec, zeromd = nocigs(H, M) 

    return eval, evec, zeromd
end


"""
LC-AGP ground state with pivot manifolds.
All pivots are assumed same and have structure
1 <= p1 < ... pk <= m. 
See DOI: 10.1063/5.0045006.

# Arguments
- `eta`: reference AGP 
- `n`: number of pairs
- `pivot`: pivot parameter
- `manifold`: pivot manifold 
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairinglcagp_pivot(eta::AbstractVector{T},
                            n::Int, 
                            pivot::T,
                            manifold::Int,
                            hamp::Float64,
                            hamtype::Symbol=:RBCS,
                            bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    manifold < 1 && error("manifold < 1 in pivotlcagp")

    # Initialize
    m = size(eta, 1)
    d = binomial(m, manifold)
    X = zeros(Int, m)
    etamat = zeros(T, d, m)

    # Get combinations
    for i = 1:m
        X[i] = i
    end
    comb = combinations(X, manifold)

    # AGPs
    i = 0
    for c in comb
	i += 1
        etamat[i, :] = pivotgen(eta, pivot, Tuple(c))
    end	

    # LC-AGP
    eval, evec, zeromd = pairinglcagp_gen(etamat, n, hamp, hamtype, bound)
    
    return eval, evec, zeromd
end
