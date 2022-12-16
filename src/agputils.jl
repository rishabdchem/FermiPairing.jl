import LinearAlgebra.BLAS: BlasFloat
using Combinatorics 

export agpoverlap, agphamoverlap, normalizeagp!, pivotgen

#======================================
Core AGP expectation values. 
======================================#

"""
<Bra|Ket>.

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
"""
function agpoverlap(U::AbstractVector{T}, 
                    V::AbstractVector{T},
	            n::Int) where T<:BlasFloat

    # Check 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in agpoverlap")

    # Initialize
    m = size(U, 1)
    X = zeros(T, m) 

    # Primary object
    X = conj(U) .* V

    # Final
    ov = getesp(X, n) 
   
    return ov 
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
- `n`: number of pairs 
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default
"""
function agphamoverlap(U::AbstractVector{T},
                       V::AbstractVector{T},
	               n::Int,
                       hamparam::Float64,
                       hamtype::Symbol=:RBCS,
                       bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in agphamoverlap")

    # Initialize
    m = size(U, 1)

    # Hamiltonian element
    if hamtype == :RBCS
        # Reduced BCS or pairing 
 
        # Initialize 
        h = zero(T) 
        G = copy(hamparam)

        # First term
        for p = 1:m
            h += ( p - G/2 ) * agptz11(U, V, n, p)
        end

        # Second term 
        for p = 1:m, q = 1:m
            p == q && continue
            h -= G * agptz02(U, V, n, (p, q,))
        end

    else
        # 1D XXZ with nearest neighbour 
    
        # Initialize 
        h = zero(T) 
        D = copy(hamparam)

        # First term
        for p = 1:m-1
            q = p+1
            h += agphpqforxxz(U, V, n, D, (p, q,)) 
        end

        # Second term 
        if bound == :PBC && m > 2 
            h += agphpqforxxz(U, V, n, D, (m, 1,)) 
        end
       
    end

    return h 
end


"""
<Bra| Hpq |Ket>, where 
Hpq = (1/2) ( P!p Pq + P!q Pp ) 
    + (D/4) ( Np Nq - Np - Nq + I ) and p =/= q. 

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
- `n`: number of pairs 
- `D`: Hamiltonian parameter 
- `hind`: Hamiltonian indices 
"""
function agphpqforxxz(U::AbstractVector{T},
                      V::AbstractVector{T},
	              n::Int,
                      D::Float64,
                      hind::Tuple) where T<:BlasFloat

    p = hind[1]
    q = hind[2]

    hpq = zero(T) 
    hpq += agptz02(U, V, n, (p, q,)) / 2 
    hpq += agptz02(U, V, n, (q, p,)) / 2
    hpq += (D / 4) * agptz22(U, V, n, (p, q,)) 
    hpq -= (D / 4) * agptz11(U, V, n, p) 
    hpq -= (D / 4) * agptz11(U, V, n, q) 
    hpq += (D / 4) * agpoverlap(U, V, n)  

    return hpq
end


"""
<Bra| P!p1 ... P!pk Pq1 ... Pqk |Ket>.
See DOI: 10.1063/1.5127850. 

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `braind`: P! indices 
- `ketind`: P indices 
"""
function agprdm(U::AbstractVector{T}, 
                V::AbstractVector{T},
	        n::Int,
	        braind::Tuple,
	        ketind::Tuple) where T<:BlasFloat

    # Check 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in agprdm")

    # Return zero if
    size(braind, 1) == size(ketind, 1) || return zero(T) 
    lenind = size(braind, 1)
    size(braind, 1) > n && return zero(T) 
    for i = 1:lenind, j = 1:i-1
        braind[i] == braind[j] && return zero(T)
        ketind[i] == ketind[j] && return zero(T)
    end

    # Initialize
    m = size(U, 1)
    X = zeros(T, m) 

    # Primary object
    X = conj(U) .* V
   
    # Prefactor 
    prefac = one(T)
    for i = 1:lenind
        i1 = braind[i]
        i2 = ketind[i]
        prefac *= conj(U[i1]) * V[i2]
        X[i1] = zero(T)
        X[i2] = zero(T)
    end

    # Final
    k = n - lenind
    rdmval = prefac * getesp(X, k) 

    return rdmval 
end


"""
Brute force computation of AGP overlaps.

# Arguments
- `etabra`: AGP bra
- `etaket`: AGP ket 
- `m`: number of levels 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function bfagpoverlaps(etabra::AbstractVector{T},
                       etaket::AbstractVector{T},
                       n::Int, 
                       hamp::Float64,
                       hamtype::Symbol=:RBCS,
                       bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    size(etabra, 1) == size(etaket, 1) || error("Fock space mismatch in bfagpoverlaps")

    # Initialize
    m = size(etabra, 1)
    d = binomial(m, n)
    X = zeros(Int, m)
    U = ones(T, d)
    V = ones(T, d)

    # Get combinations
    for j = 1:m
        X[j] = j
    end
    comb = combinations(X, n)

    # DOCI coefficients 
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            U[mu] *= etabra[pind[j]]
            V[mu] *= etaket[pind[j]]
        end    
    end 

    # Scalars 
    ov = docioverlap(U, V, m, n)
    en = docihamoverlap(U, V, m, n, hamp, hamtype, bound) 

    return en, ov 
end


#=======================================
Useful AGP expectation values. 
Np = 2 P!p Pp. 
Np^2 = 2 Np. 
=======================================#

"""
<Bra| Np |Ket>.

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `nind`: operator index
"""
function agptz11(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         nind::Int) where T<:BlasFloat

    return 2agprdm(U, V, n, (nind,), (nind,))    
end


"""
<Bra| Np Nq |Ket>.

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `nind`: operator indices 
"""
function agptz22(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         nind::Tuple) where T<:BlasFloat

    p, q = nind[1], nind[2]
    p == q && return 4agprdm(U, V, n, (p,), (p,))    

    return 4agprdm(U, V, n, (p, q,), (p, q,))    
end


"""
<Bra| Np Nq Nr |Ket>.

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `nind`: operator indices 
"""
function agptz33(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         nind::Tuple) where T<:BlasFloat

    p, q, r = nind[1], nind[2], nind[3]
    if q == r 
        return 2agptz22(U, V, n, (p, q,))    
    elseif p == r 
        return 2agptz22(U, V, n, (p, q,))    
    elseif p == q
        return 2agptz22(U, V, n, (p, r,))    
    else
        8agprdm(U, V, n, (p, q, r,), (p, q, r,))    
    end
end


"""
<Bra| P!p Pq |Ket> (p =/= q).

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `opind`: operator indices 
"""
function agptz02(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         opind::Tuple) where T<:BlasFloat

    p, q = opind[1], opind[2]
    p == q && return zero(T)

    return agprdm(U, V, n, (p,), (q,))    
end


"""
<Bra| P!p Nq Pr |Ket> (p =/= r).

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `opind`: operator indices 
"""
function agptz13(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         opind::Tuple) where T<:BlasFloat

    p, q, r = opind[1], opind[2], opind[3]
    q == r && return zero(T)
    p == r && return zero(T)
    p == q && return zero(T)

    return 2agprdm(U, V, n, (p, q,), (q, r,))    
end


"""
<Bra| P!p Nq Nr Ps |Ket> (p =/= s).

# Arguments
- `U`: bra vector
- `V`: ket vector
- `n`: number of pairs 
- `opind`: operator indices 
"""
function agptz24(U::AbstractVector{T}, 
                 V::AbstractVector{T},
	         n::Int,
	         opind::Tuple) where T<:BlasFloat

    p, q, r, s = opind[1], opind[2], opind[3], opind[4]
    p == q && return zero(T)
    p == r && return zero(T)
    p == s && return zero(T)
    q == s && return zero(T)
    r == s && return zero(T)
    
    if q == r
        return 2agptz13(U, V, n, (p, q, s,))
    else
        return 4agprdm(U, V, n, (p, q, r,), (q, r, s,))
    end
end

#======================================
Miscellaneous AGP functions.
======================================#

"""
Make <AGP|AGP> = 1.

# Arguments
- `eta`: geminal vector 
- `n`: number of pairs 
"""
function normalizeagp!(eta::AbstractVector{T}, 
                       n::Int) where T<: BlasFloat
  
    ov = agpoverlap(eta, eta, n)
    eta = eta ./ (ov ^ (1 / 2n))  

    return eta
end


"""
Update one or more AGP coefficient(s).
All pivots are assumed same.

# Arguments
- `eta`: geminal vector 
- `pivot`: pivot parameters 
- `pivind`: pivot indices 
"""
function pivotgen(eta::AbstractVector{T}, 
                  pivot::T,
                  pivind::Tuple) where T<: BlasFloat
    
    V = copy(eta)
    for i = 1:size(pivind, 1)
        V[pivind[i]] *= pivot
    end

    return V 
end


