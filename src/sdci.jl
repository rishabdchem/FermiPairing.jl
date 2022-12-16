using LinearAlgebra 
using Combinatorics 

export pairingfci, pairingcid, pairinghf 

#=====================================================
CI methods based on Slater determinants.
=====================================================#

"""
FCI.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default
- `bound`: boundary condition with PBC as default 
- `state`: state number with GS as default  
"""
function pairingfci(m::Int,
                    n::Int, 
                    hamparam::Float64,
                    hamtype::Symbol=:RBCS,
                    bound::Symbol=:PBC,
                    state::Int=1) 

    # Check 
    d = binomial(m, n)
    state > d && error("state outside FCI dimension")
 
    # Initialize
    V = zeros(Int, m, d)
    H = zeros(Float64, d, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Slater determinants
    i = 0
    for c in comb
	i += 1
	for j = 1:n    
            V[c[j], i] = 1
        end    
    end	    

    # Hamiltonian matrix
    for i = 1:d
	H[i, i] = sdhamoverlap(V[:, i], V[:, i], hamparam, hamtype, bound)
    end	    
    for i = 1:d-1, j = i+1:d
        H[i, j] = sdhamoverlap(V[:, i], V[:, j], hamparam, hamtype, bound)
        H[j, i] = H[i, j] 
    end	    

    # Diagonalize 
    @time eigval, eigvec = eigen(H)

    return eigval[state], eigvec[:, state]
end


"""
CID ground state.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default
- `bound`: boundary condition with PBC as default 
- `state`: state number with GS as default
"""
function pairingcid(m::Int,
                    n::Int, 
                    hamparam::Float64,
                    hamtype::Symbol=:RBCS,
                    bound::Symbol=:PBC,
                    state::Int=1) 

    # Check 
    d = 1 + ( n * (m - n) ) 
    state > d && error("state outside CID dimension")

    # Initialize
    X = zeros(Int, m)
    V = zeros(Int, m, d)
    H = zeros(Float64, d, d)

    # HF
    V[1:n, 1] = ones(Int, n) 

    # Unique doubles 
    p = 1
    for i = 1:n, a = n+1:m
        p += 1
        V[:, p] = phgen(m, n, (a,), (i,))
    end

    # Hamiltonian matrix
    for i = 1:d
	H[i, i] = sdhamoverlap(V[:, i], V[:, i], hamparam, hamtype, bound)
    end	    
    for i = 1:d-1, j = i+1:d
        H[i, j] = sdhamoverlap(V[:, i], V[:, j], hamparam, hamtype, bound)
        H[j, i] = H[i, j] 
    end	    

    # Solve only for ground state 
    @time eigval, eigvec = eigen(H)

    return eigval[state], eigvec[:, state]
end


"""
<HF| H |HF>, where 
|HF> = P!1 ... P!n |vac>. 

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default
- `bound`: boundary condition with PBC as default 
"""
function pairinghf(m::Int,
                   n::Int,
                   hamparam::Float64,
                   hamtype::Symbol=:RBCS,
                   bound::Symbol=:PBC) 

    V = zeros(Int, m)
    V[1:n] = ones(Int, n)

    return sdhamoverlap(V, V, hamparam, hamtype, bound) 
end

