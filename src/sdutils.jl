using Combinatorics
export sdoverlap, sdhamoverlap

#======================================
Slater determinants.
======================================#

"""
<Bra| Np |Ket>.

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
- `nind`: index p 
"""
function sdbranumket(U::AbstractVector{T}, 
                     V::AbstractVector{T},
	             nind::Int) where T<:Union{Int, Bool}

    # Check dimensions 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in sdbranumket")

    # Check if Np |Ket> = 0 
    U[nind] == 0 && return 0.0 

    # Check if <Bra| Np = 0 
    V[nind] == 0 && return 0.0 

    # Check vectors 
    U == V || return 0.0
   
    return 2.0 
end


"""
<Bra| Np Nq |Ket>.

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
- `nind`: indices
"""
function sdbranumnumket(U::AbstractVector{T}, 
                        V::AbstractVector{T},
	                nind::Tuple) where T<:Union{Int, Bool}

    # Check dimensions 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in sdbranumnumket")

    # Np Np = 2Np
    nind[1] == nind[2] && return 2sdbranumket(U, V, nind[1]) 

    # Check if Np/Nq |Ket> = 0 
    U[nind[1]] == 0 && return 0.0 
    U[nind[2]] == 0 && return 0.0 

    # Check if <Bra| Np/Nq = 0 
    V[nind[1]] == 0 && return 0.0 
    V[nind[2]] == 0 && return 0.0 

    # Check vectors 
    U == V || return 0.0
   
    return 4.0 
end


"""
<Bra| P!p Pq |Ket> (p =/= q).

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
- `pind`: indices 
"""
function sdbrapdagpket(U::AbstractVector{T}, 
                       V::AbstractVector{T},
	               pind::Tuple) where T<:Union{Int, Bool}

    # Check dimensions 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in sdbrapdagpket")

    # Check indices
    pind[1] == pind[2] && error("p = q in sdbrapdagpket")

    # Check if <Bra| P!p or Pq = 0  
    U[pind[1]] == 0 && return 0.0 
    U[pind[2]] == 1 && return 0.0 

    # Check if P!p or Pq |Ket> = 0 
    V[pind[1]] == 1 && return 0.0 
    V[pind[2]] == 0 && return 0.0 

    # Check vectors 
    X = copy(V)
    X[pind[2]] = 0 
    X[pind[1]] = 1 
    U == X || return 0.0
 
    return 1.0 
end


"""
<Bra|Ket> overlap.

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
"""
function sdoverlap(U::AbstractVector{T}, 
                   V::AbstractVector{T}) where T<:Union{Int, Bool}

    # Check dimensions 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in sdoverlap")

    # Check vectors 
    U == V || return 0.0
    
    return 1.0 
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
- `U`: bra binary vector
- `V`: ket binary vector
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function sdhamoverlap(U::AbstractVector{T},
                      V::AbstractVector{T},
                      hamparam::Float64,
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC) where T<:Union{Int, Bool}

    # Check dimensions 
    size(U, 1) == size(V, 1) || error("Fock spaces do not match in sdhamoverlap")

    # Initialize
    m = size(U, 1)

    # Hamiltonian element
    if hamtype == :RBCS
        # Reduced BCS or pairing 
   
        # Initialize 
        h = 0.0
        G = copy(hamparam)

        # First term
        for p = 1:m
            h += ( p - (G / 2) ) * sdbranumket(U, V, p)
        end

        # Second term 
        for p = 1:m, q = 1:m
            p == q && continue
            h -= G * sdbrapdagpket(U, V, (p, q,))
        end

    else
        # 1D XXZ with nearest neighbour 
    
        # Initialize 
        h = 0.0
        D = copy(hamparam)

        # First term
        for p = 1:m-1
            q = p+1
            h += sdhpqforxxz(U, V, D, (p, q,)) 
        end

        # Second term 
        if bound == :PBC && m > 2 
            h += sdhpqforxxz(U, V, D, (m, 1,)) 
        end
       
    end

    return h 
end


"""
<Bra| Hpq |Ket>, where 
Hpq = (1/2) ( P!p Pq + P!q Pp ) + (D/4) ( Np Nq - Np - Nq + I ). 

# Arguments
- `U`: bra binary vector
- `V`: ket binary vector
- `D`: Hamiltonian parameter 
- `hind`: Hamiltonian indices 
"""
function sdhpqforxxz(U::AbstractVector{T},
                     V::AbstractVector{T},
                     D::Float64,
                     hind::Tuple) where T<:Union{Int, Bool}

    p = hind[1]
    q = hind[2]

    hpq = 0.0
    hpq += sdbrapdagpket(U, V, (p, q,)) / 2 
    hpq += sdbrapdagpket(U, V, (q, p,)) / 2
    hpq += (D / 4) * sdbranumnumket(U, V, (p, q,)) 
    hpq -= (D / 4) * sdbranumket(U, V, p) 
    hpq -= (D / 4) * sdbranumket(U, V, q) 
    hpq += (D / 4) * sdoverlap(U, V)  

    return hpq
end
 
#======================================
Useful SD functions.
======================================#

"""
|SD> <- P!a1 P!a2 ... Pi1 Pi2 |HF> 

# Arguments
- `m`: number of levels 
- `n`: number of pairs 
- `virind`: P!a indices 
- `occind`: Pi indices 
"""
function phgen(m::Int,
               n::Int,
               virind::Tuple,
               occind::Tuple)

    # Check 
    size(occind, 1) == size(virind, 1) || error("Number-broken string in phgen")
    n < size(occind, 1) && error("Occupied indices more than pairs in phgen")
    m - n < size(virind, 1) && error("Virtual indices more than holes in phgen")

    # Initialize
    V = zeros(Int, m)
    V[1: n] = ones(Int, n)

    # Excitations
    for i = 1:size(occind, 1)
        V[occind[i]] = 0  
        V[virind[i]] = 1  
    end

    return V 
end


"""
|SDout> <- P!p1 ... P!pk Pq1 ... Pqk |SDin>. 

Return 0 if the operator kills |SDin>. 
Else, return DOCI index of SDout.

# Arguments
- `U`: Slater determinant vector 
- `n`: number of pairs 
- `pcind`: P! indices
- `paind`: P indices 
"""
function sdgen(U::AbstractVector{T}, 
               n::Int,
               pcind::Tuple,
               paind::Tuple) where T<:Union{Int, Bool}

    # Initialize
    m = size(U, 1) 
    k = size(paind, 1)
    V = copy(U) 

    # Check
    k == size(pcind, 1) || error("Number-breaking in sdgen")
    k > n && error("Wrong excitation rank in sdgen")
   
    # Return zero if 
    for p = 1:k
        U[paind[p]] == 0 && return 0 
    end
    for p = 1:k, q = p+1:k
        paind[p] == paind[q] && return 0 
        pcind[p] == pcind[q] && return 0 
    end

    # Pair annihilations
    for p = 1:k
        V[paind[p]] = 0
    end 

    # Pair creations 
    for p = 1:k
        V[pcind[p]] == 1 && return 0 
        V[pcind[p]] = 1
    end

    return sdindex(V, n) 
end


"""
Index in DOCI space <= Slater determinant. 
Indexing consistent with Combinatorics.jl.

# Arguments
- `V`: Slater determinant vector 
- `n`: number of pairs 
"""
function sdindex(V::AbstractVector{T}, 
                 n::Int) where T<:Union{Int, Bool}

    # Check
    m = size(V, 1)
    i = 0
    for j = 1:m 
        if V[j] == 1 
            i += 1
        end
    end
    i == n || error("Wrong number of pairs in sdindex")

    # Get combinations
    X = Vector(1:m) 
    comb = combinations(X, n)

    # Get index 
    j = 0
    for c in comb
	j += 1
        V == sdocctofull(c, m, n) && return j
    end	    

end


"""
Slater determinant <= Index in DOCI space.
Indexing consistent with Combinatorics.jl.

# Arguments
- `sdind`: Slater determinant index 
- `m`: number of levels 
- `n`: number of pairs 
"""
function invsdindex(sdind::Int,
                    m::Int,
                    n::Int) 

    # Check
    1 <= sdind <= binomial(m, n) || error("Wrong index in invsdindex")

    # Get combinations
    X = Vector(1:m) 
    comb = combinations(X, n)

    # Get index 
    V = zeros(Int, m)
    j = 0
    for c in comb
	j += 1
        j == sdind || continue
        return sdocctofull(c, m, n) 
    end	    

end


"""
Occupied orbital to full occupation representation.

# Arguments
- `U`: vector containing occ indices 
- `m`: number of levels 
- `n`: number of pairs 
"""
function sdocctofull(U::AbstractVector{T}, 
                     m::Int,
                     n::Int) where T<:Union{Int, Bool}

    # Check
    size(U, 1) == n || error("Wrong number of pairs in sdocctofull")

    # Update
    V = zeros(T, m)
    for p = 1:n
        V[U[p]] = 1
    end

    return V
end
