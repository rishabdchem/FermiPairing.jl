using FermiPairing 
using LinearAlgebra 
using Optim
using Optim: converged, iterations, minimum, minimizer

import LinearAlgebra.BLAS: BlasFloat

export pairingresagp, pairingswpconflcagp, pagpfromj1ci

#============================================
Optimizing LC-AGP.
============================================#

"""
Resonating AGP.

# Arguments
- `eta`: guess AGP vectors with m levels and n pairs 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairingresagp(eta::AbstractMatrix{T},
                       n::Int, 
                       hamp::Float64,
                       niters::Int = 1000, 
                       hamtype::Symbol=:RBCS,
                       bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)
    etalocal = lcagpmattovec(eta) 

    # Function
    function f(xvec)

        ov = lcagpoverlap(lcagpvectomat(xvec, m, d), n) 
        en = lcagphamoverlap(lcagpvectomat(xvec, m, d), n, hamp, hamtype, bound) 
        # println("elcagp inside:", ' ', en / ov)   

        return en / ov
    end
 
    # Optimize 
    res = optimize(f, etalocal, ConjugateGradient(), 
                   Optim.Options(f_tol = 1e-12,
                                 iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    # Print
    println("------------Res-AGP details-------------")
    display(res)
    println()
 
    return minimum(res), lcagpvectomat(minimizer(res), m, d) 
end


"""
LC-AGP configuration sweep.

# Arguments
- `eta`: guess AGP vectors with m levels and n pairs 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairingswpconflcagp(eta::AbstractMatrix{T},
                             n::Int, 
                             hamp::Float64,
                             swpiter::Int = 1000, 
                             optiter::Int = 5000, 
                             hamtype::Symbol=:RBCS,
                             bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)
    eps = 1e-08 
    en = randn(T)
    e0 = randn(T)
    etamat = zeros(T, d, m) 
    etamat = copy(eta) 

    # Start the loop
    for j = 1:swpiter

        # Sweep
        for a = 1:d
            en, etavec = optonelcagp(etamat, n, a, hamp, optiter, hamtype, bound) 
            etamat[a, :] = copy(etavec)
        end

        # # Print
        # println("--------")
        # println("niter: $j, en: $en")

        # Check convergence
        abs(en - e0) < eps && return en, etamat 
        e0 = copy(en)

    end
 
    # Warning
    println("Warning: did not converge in $swpiter AGP sweep iterations")

    return en, etamat
end


"""
Optimize one AGP of the LC-AGP. 

# Arguments
- `eta`: guess AGP vectors with m levels and n pairs 
- `n`: number of pairs
- `agpind`: AGP index 
- `hamp`: Hamiltonian parameter 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function optonelcagp(eta::AbstractMatrix{T},
                     n::Int, 
                     agpind::Int, 
                     hamp::Float64,
                     niters::Int = 1000, 
                     hamtype::Symbol=:RBCS,
                     bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    1 <= agpind <= size(eta, 1) || error("Wrong AGP index in optonelcagp") 
    
    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)
    etavec = zeros(T, m) 
    etamat = zeros(T, d-1, m) 
    etavec = copy(eta[agpind, :]) 
    etamat = copy(eta[[1:agpind-1; agpind+1:d], :]) 

    # Function
    function f(xvec)

        ov = lcagpoverlap(vcat(etamat[1:d-1, :], transpose(xvec)), n) 
        en = lcagphamoverlap(vcat(etamat[1:d-1, :], transpose(xvec)), n, hamp, hamtype, bound) 
        println("elcagp inside:", ' ', en / ov)   

        return en / ov
    end
 
    # Optimize 
    res = optimize(f, etavec, ConjugateGradient(), 
                   Optim.Options(f_tol = 1e-12,
                                 iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    return minimum(res), minimizer(res) 
end

#====================================
# Relevant tools.
====================================#

"""
<Psi|Psi> for a LC-AGP.

# Arguments
- `eta`: AGP vectors of m levels and n pairs
- `n`: number of pairs
- `cvec`: NOCI coefficients 
"""
function lcagpoverlap(eta::AbstractMatrix{T},
                      n::Int,
                      cvec::AbstractVector{T}=ones(T, size(eta, 1))) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)

    # Scalar
    ov = zero(T)
    for i = 1:d
	ov += ( cvec[i]^2 * agpoverlap(eta[i, :], eta[i, :], n) )
    end	    
    for i = 1:d-1, j = i+1:d
	ov += ( 2 * conj(cvec[i]) * cvec[j] * agpoverlap(eta[i, :], eta[j, :], n) )
    end

    return ov 
end


"""
<Psi| H |Psi> for a LC-AGP.

# Arguments
- `eta`: AGP vectors of m levels and n pairs
- `n`: number of pairs
- `cvec`: NOCI coefficients 
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function lcagphamoverlap(eta::AbstractArray{T},
                         n::Int,
                         hamp::Float64,
                         hamtype::Symbol=:RBCS,
                         bound::Symbol=:PBC,
                         cvec::AbstractVector{T}=ones(T, size(eta, 1))) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)
    
    # Scalar
    ov = zero(T)
    for i = 1:d
	ov += ( cvec[i]^2 * agphamoverlap(eta[i, :], eta[i, :], n, hamp, hamtype, bound) )
    end	    
    for i = 1:d-1, j = i+1:d
	ov += ( 2 * conj(cvec[i]) * cvec[j] * agphamoverlap(eta[i, :], eta[j, :], n, hamp, hamtype, bound) )
    end

    return ov 
end


"""
Pack AGP matrix into a vector.

# Argument
- `eta`: AGP vectors of m levels 
"""
function lcagpmattovec(eta::AbstractMatrix{T}) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    m = size(eta, 2)
    V = zeros(d * m) 

    # Transform
    a = 1
    b = 0
    for j = 1:d
        b = a + m - 1 
        V[a:b] = copy(eta[j, :]) 
        a = b+1
    end

    return V
end


"""
Unpack AGP coefficient vector to AGP matrix. 

# Arguments
- `eta`: LC-BTS vector 
- `m`: Number of levels 
- `d`: Number of BT states 
"""
function lcagpvectomat(eta::AbstractVector{T},
                       m::Int,
                       d::Int) where T<:BlasFloat

    # Check 
    d*m == size(eta, 1) || error("Wrong dimension in lcagpvectomat")  
    M = zeros(d, m) 

    # Transform
    a = 1
    b = 0
    for j = 1:d
        b = a + m - 1 
        M[j, :] = copy(eta[a:b])
        a = b+1
    end

    return M 
end


"""
Guess AGPs from J1-CI-AGP. 

# Arguments
- `eta`: AGP for m levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: Boundary condition with PBC as default 
- `agpnum': Number of AGPs to get  
"""
function pagpfromj1ci(eta::AbstractVector{T},
                      n::Int, 
                      hamp::Float64,
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC,
                      agpnum::Int=1) where T<:BlasFloat

    # Check 
    1 <= agpnum < size(eta, 1) || error("Wrong number of AGPs in pagpsfromj1ci") 

    # Initialize
    m = size(eta, 1)
    etamat = zeros(T, agpnum, m) 

    # Get matrices 
    H2, M2 = j1pagpoverlaps(normalizeagp!(copy(eta), n), n, hamp, hamtype, bound)
 
    # NOCI 
    @time evals, evecs, zeromd = nociall(H2, M2, agpnum+1) 

    # Get AGPs
    for a = 1:agpnum
        etamat[a, :] = expj1agp(evecs[:, a+1], eta) 
    end

    return etamat 
end

