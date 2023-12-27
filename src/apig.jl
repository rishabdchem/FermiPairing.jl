import LinearAlgebra.BLAS: BlasFloat
using Combinatorics 
using Optim 
using Optim: converged, iterations, minimum, minimizer
using Optim: only_fg!, only_fgh!

export pairingapig, papigfrompagp

#==================================================
APIG solver using DOCI subroutines.
===================================================#

"""
eagp, eta <- optimize trial APIG ground state energy.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `eta`: guess APIG 
- `niters`: number of iterations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairingapig(m::Int,
                     n::Int, 
                     hamp::Float64,
                     eta::AbstractMatrix{T},
                     niters::Int = 1000, 
                     hamtype::Symbol=:RBCS,
                     bound::Symbol=:PBC) where T<:BlasFloat

    # Function
    function f(xvec)

        en, ov = apigoverlaps(apigvectomat(xvec, m, n), hamp, hamtype, bound) 
        # println("eapig inside:", ' ', en / ov)   
 
        return en / ov
    end

    # Print
    println("----------------")
    println("Initial APIG:")
    display(eta)
    println()
    
    # Pack geminal coefficients 
    etapack = apigmattovec(eta) 

    # Optimize 
    @time res = optimize(f, etapack, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    # Print
    println("------------APIG details-------------")
    display(res)
    println()

    return minimum(res), apigvectomat(minimizer(res), m, n)
end


"""
 APIG guess from AGP.

# Arguments
- `eta`: AGP of m levels 
- `n`: number of pairs
"""
function papigfrompagp(eta::AbstractVector{T},
                       n::Int,
                       sympar::Float64=1e-03) where T<:BlasFloat
  
    # Initialize 
    m = size(eta, 1)
    G = zeros(T, n, m)

    # Guess
    for j = 1:n
        G[j, :] .= copy(eta) 
        j == 1 && continue
        for p = 1:m 
            G[j, p] += (sympar * randn(T)) 
        end
    end 

    return G 
end


"""
Brute force computation of APIG overlaps.

# Arguments
- `eta`: APIG geminal mtrix 
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function apigoverlaps(eta::AbstractMatrix{T},
                      hamp::Float64,
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    n = size(eta, 1)
    m = size(eta, 2)
    d = binomial(m, n)
    V = zeros(T, d)
    S = zeros(T, n, n)

    # DOCI coefficients 
    comb = combinations(Vector(1:m), n)
    mu = 0
    for p in comb
        mu += 1
        for j = 1:n
            S[j, :] .= copy(eta[:, p[j]]) 
        end
        V[mu] = getpermanent(S) 
    end 

    # Scalars 
    ov = docioverlap(V, V, m, n)
    en = docihamoverlap(V, V, m, n, hamp, hamtype, bound) 

    return en, ov 
end


"""
Pack APIG geminal matrix.

# Arguments
- `eta`: geminal matrix of dimensions n x m 
"""
function apigmattovec(eta::AbstractMatrix{T}) where T<: BlasFloat
  
    # Initialize 
    n = size(eta, 1)
    m = size(eta, 2)
    V = zeros(T, n * m)

    # Pack
    a = 0
    for j = 1:n, p = 1:m
        a += 1 
        V[a] = copy(eta[j, p])   
    end 

    return V 
end


"""
Unpack APIG geminal vector.

# Arguments
- `eta`: geminal vector  
- `m`: number of levels 
- `n`: number of pairs
"""
function apigvectomat(eta::AbstractVector{T},
                      m::Int,
                      n::Int) where T<: BlasFloat
  
    # Check 
    size(eta, 1) == n * m || error("Dimension mismatch in apigvectomat")
 
    # Unpack 
    M = zeros(T, n, m)
    a = 0
    for j = 1:n, p = 1:m
        a += 1 
        M[j, p] = copy(eta[a]) 
    end 

    return M 
end


