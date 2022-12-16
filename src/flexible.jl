using Combinatorics

import LinearAlgebra.BLAS: BlasFloat

export docioverlap, docihamoverlap, actopdoci 

#====================================================
Work with states that cover DOCI 

|Psi> = sum(p1 < ... < pn) C(p1 ... pn)
        P!p1 ... P!pn |->. 
====================================================#

"""
<Bra|Ket>.

# Arguments
- `U`: bra vector 
- `V`: ket vector 
- `m`: number of levels 
- `n`: number of pairs 
"""
function docioverlap(U::AbstractVector{T}, 
                     V::AbstractVector{T},
                     m::Int,
                     n::Int) where T<:BlasFloat 

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in docioverlap")
    size(U, 1) == binomial(m, n) || error("Wrong vector dimension in docioverlap")

    # Compute
    overlap = sum(conj(U) .* V) 

    return overlap
end


"""
<Bra| Np |Ket>.

# Arguments
- `U`: bra vector 
- `V`: ket vector 
- `m`: number of levels 
- `n`: number of pairs 
- `nind`: index p 
"""
function docibranumket(U::AbstractVector{T}, 
                       V::AbstractVector{T},
                       m::Int,
                       n::Int,
                       nind::Int) where T<:BlasFloat 

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in docibranumket")
    size(U, 1) == binomial(m, n) || error("Wrong vector dimension in docibranumket")

    # Initialize
    d = binomial(m, n)
    S = zeros(Int, m, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Slater determinants
    i = 0
    for c in comb
	i += 1
	for j = 1:n    
            S[c[j], i] = 1
        end    
    end	    

    # Compute
    overlap = 0.0 
    for i = 1:d, j = 1:d
        overlap += conj(U[i]) * V[j] * sdbranumket(S[:, i], S[:, j], nind)    
    end

    return overlap
end


"""
<Bra| Np Nq |Ket>.

# Arguments
- `U`: bra vector 
- `V`: ket vector 
- `m`: number of levels 
- `n`: number of pairs 
- `nind`: indices 
"""
function docibranumnumket(U::AbstractVector{T}, 
                          V::AbstractVector{T},
                          m::Int,
                          n::Int,
                          nind::Tuple) where T<:BlasFloat 

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in docibranumnumket")
    size(U, 1) == binomial(m, n) || error("Wrong vector dimension in docibranumnumket")
    nind[1] == nind[2] && return 2docibranumket(U, V, m, n, nind[1]) 

    # Initialize
    d = binomial(m, n)
    S = zeros(Int, m, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Slater determinants
    i = 0
    for c in comb
	i += 1
	for j = 1:n    
            S[c[j], i] = 1
        end    
    end	    

    # Compute
    overlap = 0.0 
    for i = 1:d, j = 1:d
        overlap += conj(U[i]) * V[j] * sdbranumnumket(S[:, i], S[:, j], nind)    
    end

    return overlap
end


"""
<Bra| P!p Pq |Ket> (p =/= q).

# Arguments
- `U`: bra vector 
- `V`: ket vector 
- `m`: number of levels 
- `n`: number of pairs 
- `pind`: indices 
"""
function docibrapdagpket(U::AbstractVector{T}, 
                         V::AbstractVector{T},
                         m::Int,
                         n::Int,
                         pind::Tuple) where T<:BlasFloat 

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in docibrapdagpket")
    size(U, 1) == binomial(m, n) || error("Wrong vector dimension in docibrapdagpket")
    pind[1] == pind[2] && return docibranumket(U, V, m, n, pind[1]) / 2 

    # Initialize
    d = binomial(m, n)
    S = zeros(Int, m, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Slater determinants
    i = 0
    for c in comb
	i += 1
	for j = 1:n    
            S[c[j], i] = 1
        end    
    end	    

    # Compute
    overlap = 0.0 
    for i = 1:d, j = 1:d
        overlap += conj(U[i]) * V[j] * sdbrapdagpket(S[:, i], S[:, j], pind)    
    end

    return overlap
end


"""
<Bra| H |Ket>.

RBCS: 
H = sum(p) ( p - G/2 ) Np 
  - G sum(p =/= q) P!p Pq.

1D XXZ:  
H = specialsum(pq) (1/2) ( P!p Pq + P!q Pp ) 
  + (D/4) ( Np Nq - Np - Nq + I ), 
where "specialsum" indicates NN with boundary condition.

# Arguments
- `U`: bra vector 
- `V`: ket vector 
- `m`: number of levels 
- `n`: number of pairs 
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function docihamoverlap(U::AbstractVector{T}, 
                        V::AbstractVector{T},
                        m::Int,
                        n::Int,
                        hamparam::Float64,
                        hamtype::Symbol=:RBCS,
                        bound::Symbol=:PBC) where T<:BlasFloat 

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in docihamoverlap")
    size(U, 1) == binomial(m, n) || error("Wrong vector dimension in docihamoverlap")

    # Initialize
    d = binomial(m, n)
    S = zeros(Int, m, d)
    H = zeros(Float64, d, d)

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Slater determinants
    i = 0
    for c in comb
	i += 1
	for j = 1:n    
            S[c[j], i] = 1
        end    
    end	    

    # Hamiltonian matrix
    for i = 1:d
	H[i, i] = sdhamoverlap(S[:, i], S[:, i], hamparam, hamtype, bound)
    end	    
    for i = 1:d-1, j = i+1:d
        H[i, j] = sdhamoverlap(S[:, i], S[:, j], hamparam, hamtype, bound)
        H[j, i] = H[i, j] 
    end	

    # Compute
    overlap = 0.0 
    for i = 1:d, j = 1:d
        overlap += conj(U[i]) * V[j] * H[i, j] 
    end

    return overlap
end


"""
P!p1 ... P!pk Pq1 ... Pqk |Psi>. 

This algorithm gets help from the fact that the operator string
does not map many-to-one for the basis SDs.

# Arguments
- `U`: DOCI coefficients 
- `m`: number of levels 
- `n`: number of pairs 
- `pcind`: P! indices 
- `paind`: P indices 
"""
function actopdoci(U::AbstractVector{T}, 
                   m::Int,
                   n::Int,
                   pcind::Tuple,
                   paind::Tuple) where T<:BlasFloat 

    # Check
    d = binomial(m, n)
    size(U, 1) == d || error("Wrong vector dimension in actopdoci")
    size(paind, 1) == size(pcind, 1) || error("Number-breaking in actopdoci")
    size(paind, 1) > n && error("Wrong excitation rank in actopdoci")

    # Initialize
    X = Vector(1:m) 
    Y = zeros(Int, m) 
    V = zeros(T, d)

    # Get combinations
    comb = combinations(X, n)

    # Loop over all basis states 
    i = 0
    for c in comb
	i += 1
        # Slater determinant 
	for j = 1:n    
            Y[c[j]] = 1
        end    
        # Act with operators
        indout = sdgen(Y, n, pcind, paind) 
        # Update V 
        if indout != 0
            V[indout] = copy(U[i])  
        end
        # Reset SD 
        Y .= 0
    end	

    return V 
end


