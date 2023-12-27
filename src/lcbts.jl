using LinearAlgebra 
using Combinatorics 

import LinearAlgebra.BLAS: BlasFloat

export genlcbts, btsbraket, btsbrahamket

#=====================================
LC-BTS. 
=====================================#

"""
LC-BTS ground state with fixed BT states.

# Arguments
- `eta`: bts matrices of m levels and n pairs
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function genlcbts(eta::AbstractArray{T},
                  hamp::Float64,
                  hamtype::Symbol=:RBCS,
                  bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    M = zeros(T, d, d)
    H = zeros(T, d, d)
    etalocal = copy(eta)

    # Normalize BT states
    for i = 1:d
        etalocal[i, :, :] = normalizebts!(etalocal[i, :, :])
    end

    # Matrices 
    for i = 1:d
	M[i, i] = btsbraket(etalocal[i, :, :], etalocal[i, :, :])
	H[i, i] = btsbrahamket(etalocal[i, :, :], etalocal[i, :, :], hamp, hamtype, bound)
    end	    
    for i = 1:d-1, j = i+1:d
	M[i, j] = btsbraket(etalocal[i, :, :], etalocal[j, :, :])
        M[j, i] = conj(M[i, j]) 
        H[i, j] = btsbrahamket(etalocal[i, :, :], etalocal[j, :, :], hamp, hamtype, bound)
        H[j, i] = conj(H[i, j]) 
    end	    

    # Solve for ground state 
    @time eval, evec, zeromd = nocigs(H, M) 

    return eval, evec, zeromd
end

#=============================================
Core BTS transition overlaps.
=============================================#

"""
<Bra|Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
"""
function btsbraket(etabra::AbstractMatrix{T},
                   etaket::AbstractMatrix{T}) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsbraket")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && error("n > m in btsbraket")

    # Trivial 
    m == 0 && return one(T) 
    n == 0 && return one(T) 

    return getbtp(conj.(etabra) .* etaket) 
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
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default
"""
function btsbrahamket(etabra::AbstractMatrix{T},
                      etaket::AbstractMatrix{T},
                      hamparam::Float64,
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsbrahamket")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && error("n > m in btsbrahamket")

    # Hamiltonian element
    if hamtype == :RBCS
        # Reduced BCS or pairing 
   
        # Initialize 
        h = 0.0
        G = copy(hamparam)

        # First term
        for p = 1:m
            h += ( p - (G / 2) ) * btstz11(etabra, etaket, p)
        end

        # Second term 
        for p = 1:m, q = 1:m
            p == q && continue
            h -= G * btstz02(etabra, etaket, (p, q,))
        end

    else
        # 1D XXZ with nearest neighbour 
    
        # Initialize 
        h = 0.0
        D = copy(hamparam)

        # First term
        for p = 1:m-1
            q = p+1
            h += btsbrahpqketforxxz(etabra, etaket, D, (p, q,)) 
        end

        # Second term 
        if bound == :PBC && m > 2 
            h += btsbrahpqketforxxz(etabra, etaket, D, (m, 1,)) 
        end
       
    end

    return h 
end


"""
<Bra| Hpq |Ket>, where 
Hpq = (1/2) ( P!p Pq + P!q Pp ) + (D/4) ( Np Nq - Np - Nq + I ). 

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs
- `D`: Hamiltonian parameter 
- `hind`: Hamiltonian indices 
"""
function btsbrahpqketforxxz(etabra::AbstractMatrix{T},
                            etaket::AbstractMatrix{T},
                            D::Float64,
                            hind::Tuple) where T<:BlasFloat

    # Initialize
    p = hind[1]
    q = hind[2]
    
    # Compute
    hpq = zero(T) 
    hpq += btstz02(etabra, etaket, (p, q,)) / 2 
    hpq += btstz02(etabra, etaket, (q, p,)) / 2
    hpq += (D / 4) * btstz22(etabra, etaket, (p, q,)) 
    hpq -= (D / 4) * btstz11(etabra, etaket, p) 
    hpq -= (D / 4) * btstz11(etabra, etaket, q) 
    hpq += (D / 4) * btsbraket(etabra, etaket)

    return hpq
end


"""
Brute force computation of transition overlaps.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function bfbtstransitions(etabra::AbstractMatrix{T},
                          etaket::AbstractMatrix{T},
                          hamp::Float64,
                          hamtype::Symbol=:RBCS,
                          bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in bfbtstransitions")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && error("n > m in bfbtstransitions")

    # Initialize
    d = binomial(m, n)
    U = ones(T, d)
    V = ones(T, d)

    # Get combinations
    X = Vector(1:m) 
    comb = combinations(X, n)

    # Coefficients 
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            U[mu] *= etabra[j, pind[j]]
            V[mu] *= etaket[j, pind[j]]
        end    
    end 
 
    # Scalars 
    ov = docioverlap(U, V, m, n)
    en = docihamoverlap(U, V, m, n, hamp, hamtype, bound) 

    return en, ov 
end

#=========================================
BTS transition RDMs. 
=========================================#

"""
<Bra| Np |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `nind`: operator index 
"""
function btstz11(etabra::AbstractMatrix{T}, 
                 etaket::AbstractMatrix{T}, 
	         nind::Int) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btstz11")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && return zero(T) 
    nind > m && return zero(T) 

    # Initialize
    p = nind
    X = conj(etabra) .* etaket

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
<Bra| Np Nq |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `nind`: operator indices 
"""
function btstz22(etabra::AbstractMatrix{T}, 
                 etaket::AbstractMatrix{T}, 
	         nind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btstz22")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && error("n > m in btsnumnumbts")
    nind[1] == nind[2] && return 2btstz11(etabra, etaket, nind[1]) 

    # Initialize
    if nind[1] > nind[2]
        p = nind[1]
        q = nind[2]
    else
        p = nind[2]
        q = nind[1]
    end
    X = conj.(etabra) .* etaket

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
<Bra| P!p Pq |Ket>.
We use an iterative scheme like SumBTP.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `pind`: operator indices 
"""
function btstz02(etabra::AbstractMatrix{T}, 
                 etaket::AbstractMatrix{T}, 
	         pind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btstz02")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Trivial
    n > m && return zero(T) 
    intscal = zero(T)
    pind[1] == pind[2] && return btstz11(etabra, etaket, pind[1]) / 2
    n == 1 && return conj(etabra[1, pind[1]])*etaket[1, pind[2]]  
    if pind[1] == m 
        intscal = btsasympa(etabra[1:n-1, 1:m-1], etaket[1:n, 1:m-1], pind[2])
        return conj(etabra[n, m])*intscal
    elseif pind[2] == m 
        intscal = btsasympc(etabra[1:n, 1:m-1], etaket[1:n-1, 1:m-1], pind[1])
        return etaket[n, m]*intscal
    end

    # Initialize
    p, q = pind[1], pind[2]
    D = zeros(T, m, n) 
    D[:, 1] .= conj(etabra[1, p])*etaket[1, q]  
    if p < q 

        # Init
        for j = btprowmin(m, n, q):btprowmax(n, q) 
            j == 1 && continue
            intscal = btsasympc(etabra[1:j, 1:q-1], etaket[1:j-1, 1:q-1], p)     
            D[q, j] = ( etaket[j, q] * intscal ) 
        end 

        # Recursion 
        for r = q+1:m, j = btprowmin(m, n, r):btprowmax(n, r) 
            j == 1 && continue
            D[r, j] = D[r-1, j]
            D[r, j] += ( conj(etabra[j, r]) * etaket[j, r] * D[r-1, j-1] )
        end 

    else 

        # Init
        for j = btprowmin(m, n, p):btprowmax(n, p) 
            j == 1 && continue
            intscal = btsasympa(etabra[1:j-1, 1:p-1], etaket[1:j, 1:p-1], q)     
            D[p, j] = ( conj(etabra[j, p]) * intscal ) 
        end

        # Recursion 
        for r = p+1:m, j = btprowmin(m, n, r):btprowmax(n, r) 
            j == 1 && continue
            D[r, j] = D[r-1, j]
            D[r, j] += ( conj(etabra[j, r]) * etaket[j, r] * D[r-1, j-1] )
        end 

    end

    return D[m, n] 
end


"""
<Bra(m, n)| P!p |Ket(m, n-1)>.
We use an iterative scheme like SumESP.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n-1 pairs
- `pind`: operator index
"""
function btsasympc(etabra::AbstractMatrix{T}, 
                   etaket::AbstractMatrix{T}, 
	           pind::Int) where T<:BlasFloat

    # Check 
    size(etabra, 2) == size(etaket, 2) || error("Bra-ket dimension mismatch in btsasympc")
    size(etabra, 1) == size(etaket, 1) + 1 || error("Bra-ket dimension mismatch in btsasympc")
    m = size(etabra, 2)
    n = size(etabra, 1)

    # Trivial cases 
    n > m && return zero(T) 
    n == 1 && return conj(etabra[1, pind]) 
    intscal = zero(T)
    if pind == m 
        intscal = btsbraket(etabra[1:n-1, 1:m-1], etaket[1:n-1, 1:m-1])
        return conj(etabra[n, m])*intscal    
    end

    # Initialize
    D = zeros(T, m, n) 
    D[:, 1] .= conj(etabra[1, pind]) 
    if pind > 1 
        for j = btprowmin(m, n, pind):btprowmax(n, pind) 
            j == 1 && continue
            intscal = btsbraket(etabra[1:j-1, 1:pind-1], etaket[1:j-1, 1:pind-1]) 
            D[pind, j] = ( conj(etabra[j, pind]) * intscal )
        end 
    end 
    # Recursion 
    for p = pind+1:m, j = btprowmin(m, n, p):btprowmax(n, p) 
        j == 1 && continue
        D[p, j] = D[p-1, j]
        D[p, j] += ( conj(etabra[j, p]) * etaket[j-1, p] * D[p-1, j-1] )
    end 

    return D[m, n]
end


"""
<Bra(m, n-1)| Pp |Ket(m, n)>.
We use an iterative scheme like SumESP.

# Arguments
- `etabra`: bra BTS of m levels and n-1 pairs 
- `etaket`: ket BTS of m levels and n pairs
- `pind`: operator index
"""
function btsasympa(etabra::AbstractMatrix{T}, 
                   etaket::AbstractMatrix{T}, 
	           pind::Int) where T<:BlasFloat

    # Check 
    size(etabra, 2) == size(etaket, 2) || error("Bra-ket dimension mismatch in btsasympa")
    size(etabra, 1) == size(etaket, 1) - 1 || error("Bra-ket dimension mismatch in btsasympa")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Trivial cases 
    n > m && return zero(T) 
    n == 1 && return etaket[1, pind] 
    intscal = zero(T)
    if pind == m 
        intscal = btsbraket(etabra[1:n-1, 1:m-1], etaket[1:n-1, 1:m-1])
        return etaket[n, m]*intscal    
    end

    # Initialize
    D = zeros(T, m, n) 
    D[:, 1] .= etaket[1, pind] 
    if pind > 1 
        for j = btprowmin(m, n, pind):btprowmax(n, pind) 
            j == 1 && continue
            intscal = btsbraket(etabra[1:j-1, 1:pind-1], etaket[1:j-1, 1:pind-1]) 
            D[pind, j] = ( etaket[j, pind] * intscal )
        end 
    end 
    # Recursion 
    for p = pind+1:m, j = btprowmin(m, n, p):btprowmax(n, p) 
        j == 1 && continue
        D[p, j] = D[p-1, j]
        D[p, j] += ( conj(etabra[j-1, p]) * etaket[j, p] * D[p-1, j-1] )
    end 

    return D[m, n] 
end

#==================================================
Some special RDMs.
We use the relation 

|BTS'> = bNp |BTS>, where 

bNp = Np + beta,
eta'(:, p) = (1 + 2/beta) eta(:, p).
==================================================#

"""
<Bra| Np P!q Pr Ns |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `opind`: operator indices 
"""
function btsnumz02num(etabra::AbstractMatrix{T}, 
                      etaket::AbstractMatrix{T}, 
	              opind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsnumz02num")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Initialize
    beta = 2
    p, q, r, s = opind[1], opind[2], opind[3], opind[4]
    U = copy(etabra) 
    V = copy(etaket) 

    # <P!q Pr>
    rdm1 = btstz02(etabra, etaket, (q, r,))

    # <Np P!q Pr> and <P!q Pr Ns>
    rdm2 = btsnumz02(etabra, etaket, (p, q, r,))
    rdm3 = btsz02num(etabra, etaket, (q, r, s,))

    # <bNp P!q Pr bNs>
    U[:, p] .*= (1 + (2 / beta))  
    V[:, s] .*= (1 + (2 / beta))  
    rdmval = ( (beta ^ 2) * btstz02(U, V, (q, r,)) )

    # Final
    rdmval -= ( (beta ^ 2) * rdm1 ) 
    rdmval -= ( beta * rdm2 ) 
    rdmval -= ( beta * rdm3 ) 
    
    return rdmval 
end


"""
<Bra| Np P!q Pr |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `opind`: operator indices 
"""
function btsnumz02(etabra::AbstractMatrix{T}, 
                   etaket::AbstractMatrix{T}, 
	           opind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsnumz02num")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Initialize
    beta = 2
    p, q, r = opind[1], opind[2], opind[3]
    U = copy(etabra) 

    # <P!q Pr>
    rdmint = btstz02(etabra, etaket, (q, r,))

    # <bNp P!q Pr>
    U[:, p] .*= (1 + (2 / beta))  
    rdmval = btstz02(U, etaket, (q, r,))

    # Final
    rdmval -= rdmint 
    rdmval *= beta 
    
    return rdmval 
end


"""
<Bra| P!p Pq Nr |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `opind`: operator indices 
"""
function btsz02num(etabra::AbstractMatrix{T}, 
                   etaket::AbstractMatrix{T}, 
	           opind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsnumz02num")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Initialize
    beta = 2
    p, q, r = opind[1], opind[2], opind[3]
    V = copy(etaket) 

    # <P!p Pq>
    rdmint = btstz02(etabra, etaket, (p, q,))

    # <P!p Pq bNr>
    V[:, r] .*= (1 + (2 / beta))  
    rdmval = btstz02(etabra, V, (p, q,))

    # Final
    rdmval -= rdmint  
    rdmval *= beta 
    
    return rdmval 
end


"""
<Bra| Np Nq Nr Ns |Ket>.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `opind`: operator indices 
"""
function btstz44(etabra::AbstractMatrix{T}, 
                 etaket::AbstractMatrix{T}, 
	         opind::Tuple) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsnumz02num")
    m = size(etaket, 2)
    n = size(etaket, 1)

    # Trivial 
    p, q, r, s = opind[1], opind[2], opind[3], opind[4]
    if r == s
        return 4btsnumz02num(etabra, etaket, (p, q, q, r,))
    elseif q == s
        return 4btsnumz02num(etabra, etaket, (p, q, q, r,))
    elseif p == s
        return 4btsnumz02num(etabra, etaket, (p, q, q, r,))
    elseif q == r
        return 4btsnumz02num(etabra, etaket, (p, q, q, s,))
    elseif p == r
        return 4btsnumz02num(etabra, etaket, (p, q, q, s,))
    elseif p == q
        return 4btsnumz02num(etabra, etaket, (p, r, r, s,))
    end

    # Initialize
    beta = 2
    U = zeros(T, n, m)
    V = zeros(T, n, m)
    U .= copy(etabra) 
    V .= copy(etaket) 

    # <Nq Nr>
    rdm1 = btstz22(etabra, etaket, (q, r,))

    # <Np Nq Nr> and <Nq Nr Ns>
    rdm2 = 2btsnumz02num(etabra, etaket, (p, q, q, r,))
    rdm3 = 2btsnumz02num(etabra, etaket, (q, r, r, s,))
    
    # <bNp Nq Nr bNs>
    U[:, p] .*= (1 + (2 / beta))  
    V[:, s] .*= (1 + (2 / beta))  
    rdmval = ( (beta ^ 2) * btstz22(U, V, (q, r,)) )
    
    # Final
    rdmval -= ( (beta^2) * rdm1 ) 
    rdmval -= ( beta * rdm2 ) 
    rdmval -= ( beta * rdm3 ) 

    return rdmval 
end


"""
<Bra| Np H Nq |Ket>.

RBCS: 
H = sum(p) ( p - G/2 ) Np 
  - G sum(p =/= q) P!p Pq.

1D XXZ:  
H = specialsum(pq) (1/2) ( P!p Pq + P!q Pp ) 
  + (D/4) ( Np Nq - Np - Nq + I ), 
where "specialsum" indicates NN with boundary condition.

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs 
- `opind`: operator indices 
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default
"""
function btsnumbrahamnumket(etabra::AbstractMatrix{T},
                            etaket::AbstractMatrix{T},
                            opind::Tuple,
                            hamparam::Float64,
                            hamtype::Symbol=:RBCS,
                            bound::Symbol=:PBC) where T<:BlasFloat

    # Check 
    size(etabra) == size(etaket) || error("Bra-ket dimension mismatch in btsbrahamket")
    m = size(etaket, 2)
    n = size(etaket, 1)
    n > m && error("n > m in btsbrahamket")

    # Hamiltonian element
    p, q = opind[1], opind[2]
    if hamtype == :RBCS
        # Reduced BCS or pairing 
   
        # Initialize 
        h = 0.0
        G = copy(hamparam)

        # First term
        for p0 = 1:m
            h += ( 2p0 - G ) * btsnumz02num(etabra, etaket, (p, p0, p0, q,))
        end

        # Second term 
        for p0 = 1:m, p1 = 1:m
            p0 == p1 && continue
            h -= G * btsnumz02num(etabra, etaket, (p, p0, p1, q,))
        end

    else
        # 1D XXZ with nearest neighbour 
    
        # Initialize 
        h = 0.0
        D = copy(hamparam)

        # First term
        for p0 = 1:m-1
            p1 = p0 + 1
            h += btsbranumhpqnumketforxxz(etabra, etaket, D, (p, p0, p1, q,)) 
        end

        # Second term 
        if bound == :PBC && m > 2 
            h += btsbranumhpqnumketforxxz(etabra, etaket, D, (p, q,), (m, 1,)) 
        end
       
    end

    return h 
end


"""
<Bra| Na Hpq Nb |Ket>, where 
Hpq = (1/2) ( P!p Pq + P!q Pp ) + (D/4) ( Np Nq - Np - Nq + I ). 

# Arguments
- `etabra`: bra BTS of m levels and n pairs 
- `etaket`: ket BTS of m levels and n pairs
- `D`: Hamiltonian parameter 
- `hind`: Hamiltonian indices 
"""
function btsbranumhpqnumketforxxz(etabra::AbstractMatrix{T},
                                  etaket::AbstractMatrix{T},
                                  D::Float64,
                                  opind::Tuple,
                                  hind::Tuple) where T<:BlasFloat

    # Initialize
    a, b = opind[1], opind[2]
    p, q = hind[1], hind[2]
    
    # Compute
    hpq = zero(T) 
    hpq += btsnumz02num(etabra, etaket, (a, p, q, b,)) / 2 
    hpq += btsnumz02num(etabra, etaket, (a, q, p, b,)) / 2
    hpq += (D / 4) * btstz44(etabra, etaket, (a, p, q, b,)) 
    hpq -= (D / 2) * btsnumz02num(etabra, etaket, (a, p, p, b,)) 
    hpq -= (D / 2) * btsnumz02num(etabra, etaket, (a, q, q, b,)) 
    hpq += (D / 4) * btstz22(etabra, etaket, (a, b,))

    return hpq
end


