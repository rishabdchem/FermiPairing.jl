using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export btsoverlap, btshamoverlap, normalizebts!

#==========================================
BTS RDMs. 
==========================================#

"""
<BTS|BTS>.

# Argument
- `eta`: BTS matrix of m levels and n pairs 
"""
function btsoverlap(eta::AbstractMatrix{T}) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)
    n > m && error("n > m in btsoverlap")

    # Trivial 
    m == 0 && return one(T) 
    n == 0 && return one(T) 

    return getbtp(conj.(eta) .* eta) 
end


"""
Make <BTS|BTS> = 1.

# Argument
- `eta`: BTS matrix of m levels and n pairs 
"""
function normalizebts!(eta::AbstractMatrix{T}) where T<: BlasFloat
  
    n = size(eta, 1)
    ov = btsoverlap(eta)
    eta = eta ./ (ov ^ (1 / 2n))  

    return eta
end


"""
<BTS| H |BTS>.

RBCS: 
H = sum(p) ( eps[p] - G/2 ) Np 
  - G sum(p =/= q) P!p Pq.

1D XXZ:  
H = specialsum(pq) (1/2) ( P!p Pq + P!q Pp ) 
  + (D/4) ( Np Nq - Np - Nq + I ), 
where "specialsum" indicates NN with boundary condition.

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default
- `eps`: pairing levels 
"""
function btshamoverlap(eta::AbstractMatrix{T},
                       hamparam::Float64,
                       hamtype::Symbol=:RBCS,
                       bound::Symbol=:PBC,
                       eps::Vector{Int}=Vector(1:size(eta, 2))) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)
    n > m && error("n > m in btshamoverlap")

    # Hamiltonian element
    if hamtype == :RBCS
        # Reduced BCS or pairing 
   
        # Initialize 
        h = 0.0
        G = copy(hamparam)

        # First term
        for p = 1:m
            h += ( eps[p] - (G / 2) ) * btsz11(eta, p)
        end

        # Second term 
        for p = 1:m, q = 1:m
            p == q && continue
            h -= G * btsz02(eta, (p, q,))
        end

    else
        # 1D XXZ with nearest neighbour 
    
        # Initialize 
        h = 0.0
        D = copy(hamparam)

        # First term
        for p = 1:m-1
            q = p+1
            h += btshpqbtsforxxz(eta, D, (p, q,)) 
        end

        # Second term 
        if bound == :PBC && m > 2 
            h += btshpqbtsforxxz(eta, D, (m, 1,)) 
        end
       
    end

    return h 
end


"""
<BTS| Hpq |BTS>, where 
Hpq = (1/2) ( P!p Pq + P!q Pp ) + (D/4) ( Np Nq - Np - Nq + I ). 

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `D`: Hamiltonian parameter 
- `hind`: Hamiltonian indices 
"""
function btshpqbtsforxxz(eta::AbstractMatrix{T},
                         D::Float64,
                         hind::Tuple) where T<:BlasFloat

    # Initialize
    p = hind[1]
    q = hind[2]
    
    # Compute
    hpq = zero(T) 
    hpq += btsz02(eta, (p, q,)) / 2 
    hpq += btsz02(eta, (q, p,)) / 2
    hpq += (D / 4) * btsz22(eta, (p, q,)) 
    hpq -= (D / 4) * btsz11(eta, p) 
    hpq -= (D / 4) * btsz11(eta, q) 
    hpq += (D / 4) * btsoverlap(eta)

    return hpq
end


"""
Brute force computation of overlaps.

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function bfbtsoverlaps(eta::AbstractMatrix{T},
                       hamp::Float64,
                       hamtype::Symbol=:RBCS,
                       bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    m = size(eta, 2)
    n = size(eta, 1)
    n > m && error("n > m in bfbtsoverlaps")

    # Initialize
    d = binomial(m, n)
    V = ones(T, d)

    # Get combinations
    X = Vector(1:m) 
    comb = combinations(X, n)

    # Coefficients 
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            V[mu] *= eta[j, pind[j]]
        end    
    end 
 
    # Scalars 
    ov = docioverlap(V, V, m, n)
    en = docihamoverlap(V, V, m, n, hamp, hamtype, bound) 

    return en, ov 
end


"""
<BTS| Np |BTS>.
We use the general BTS partition. 

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `nind`: operator index 
"""
function btsz11(eta::AbstractMatrix{T}, 
	        nind::Int) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)
    n > m && return zero(T) 
    nind > m && return zero(T) 

    # Initialize
    p = nind
    X = conj.(eta) .* eta

    # <P!p Pp>    
    rdmval = zero(T)
    btpint = one(T)
    for j = btprowmin(m, n, p):btprowmax(n, p)
        # Intermediates
        btpint *= getbtp(X[1:j-1, 1:p-1])
        btpint *= getbtp(X[j+1:n, p+1:m])
        # RDM
        rdmval += ( X[j, p] * btpint ) 
        # Update
        btpint = one(T)
    end

    return 2rdmval 
end


"""
<BTS| Np Nq |BTS>.
We use the general BTS partition. 

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `nind`: operator indices 
"""
function btsz22(eta::AbstractMatrix{T}, 
	        nind::Tuple) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)
    n > m && error("n > m in btsnumnumbts")
    nind[1] == nind[2] && return 2btsz11(eta, nind[1]) 

    # Initialize
    if nind[1] > nind[2]
        p = nind[1]
        q = nind[2]
    else
        p = nind[2]
        q = nind[1]
    end
    X = conj.(eta) .* eta

    # <P!p Pp P!q Pq>    
    rdmval = zero(T)
    btpint = one(T)
    for j = btprowmin(m, n, p):btprowmax(n, p)
        for k = btprowmin(p-1, j-1, q):btprowmax(j-1, q)
            # Intermediates
            btpint *= getbtp(X[1:k-1, 1:q-1]) 
            btpint *= getbtp(X[k+1:j-1, q+1:p-1]) 
            btpint *= getbtp(X[j+1:n, p+1:m]) 
            # RDM
            rdmval += ( X[j, p] * X[k, q] * btpint ) 
            # Update
            btpint = one(T)
        end
    end

    return 4rdmval 
end


"""
<BTS| P!p Pq |BTS>.
We use an iterative scheme like SumBTP.

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `pind`: operator indices 
"""
function btsz02(eta::AbstractMatrix{T}, 
	        pind::Tuple) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)

    # Trivial
    n > m && return zero(T) 
    pind[1] == pind[2] && return btsz11(eta, pind[1]) / 2
    n == 1 && return conj(eta[1, pind[1]])*eta[1, pind[2]]  
    pind[1] == m && return conj(eta[n, m])*btsasympa(eta[1:n, 1:m-1], pind[2])
    pind[2] == m && return eta[n, m]*btsasympc(eta[1:n, 1:m-1], pind[1])

    # Initialize
    if pind[1] > pind[2]
        p = pind[2]
        q = pind[1]
    else
        p = pind[1]
        q = pind[2]
    end
    D = zeros(T, m, n) 
    D[:, 1] .= conj(eta[1, p])*eta[1, q]  
    for j = btprowmin(m, n, q):btprowmax(n, q) 
        j == 1 && continue
        D[q, j] = ( eta[j, q] * btsasympc(eta[1:j, 1:q-1], p) )    
    end 
    # Recursion 
    for r = q+1:m, j = btprowmin(m, n, r):btprowmax(n, r) 
        j == 1 && continue
        D[r, j] = D[r-1, j]
        D[r, j] += ( eta[j, r]^2 * D[r-1, j-1] )
    end 

    pind[1] > pind[2] && return conj(D[m, n])
    return D[m, n] 
end


"""
<BTS(m, n)| P!p |BTS(m, n-1)>.
We use an iterative scheme like SumESP.

# Arguments
- `eta`: BTS matrix of n levels and m pairs 
- `pind`: operator index
"""
function btsasympc(eta::AbstractMatrix{T}, 
	           pind::Int) where T<:BlasFloat

    return conj(btsasympa(eta, pind)) 
end


"""
<BTS(m, n-1)| Pp |BTS(m, n)>.
We use an iterative scheme like SumESP.

# Arguments
- `eta`: BTS matrix of m levels and n pairs 
- `pind`: operator index
"""
function btsasympa(eta::AbstractMatrix{T}, 
	           pind::Int) where T<:BlasFloat

    # Check 
    m = size(eta, 2)
    n = size(eta, 1)

    # Trivial cases 
    n > m && return zero(T) 
    n == 1 && return eta[1, pind] 
    pind == m && return eta[n, m]*btsoverlap(eta[1:n-1, 1:m-1])

    # Initialize
    D = zeros(T, m, n) 
    D[:, 1] .= eta[1, pind] 
    if pind > 1 
        for j = btprowmin(m, n, pind):btprowmax(n, pind) 
            j == 1 && continue
            D[pind, j] = ( eta[j, pind] * btsoverlap(eta[1:j-1, 1:pind-1]) )
        end 
    end 
    # Recursion 
    for p = pind+1:m, j = btprowmin(m, n, p):btprowmax(n, p) 
        j == 1 && continue
        D[p, j] = D[p-1, j]
        D[p, j] += ( conj(eta[j-1, p]) * eta[j, p] * D[p-1, j-1] )
    end 

    return D[m, n] 
end


