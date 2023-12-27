using LinearAlgebra 
using Optim
using Optim: converged, iterations, minimum, minimizer

import LinearAlgebra.BLAS: BlasFloat

export resbts, swpconflcbts 

#=====================================
Optimizing LC-BTS.
=====================================#

"""
Optimize trial LC-BTS energy.

# Arguments
- `eta`: guess bts matrices with m levels and n pairs 
- `hamp`: Hamiltonian parameter 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function resbts(eta::AbstractArray{T},
                hamp::Float64,
                niters::Int = 1000, 
                hamtype::Symbol=:RBCS,
                bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    etalocal = lcbtstentovec(eta) 

    # Function
    function f(xvec)

        ov = lcbtsoverlap(lcbtsvectoten(xvec, m, n, d)) 
        en = lcbtshamoverlap(lcbtsvectoten(xvec, m, n, d), hamp, hamtype, bound) 
        # println("elcbts inside:", ' ', en / ov)   

        return en / ov
    end
 
    # Optimize 
    @time res = optimize(f, etalocal, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    # Print
    println("------------Res-BTS details-------------")
    display(res)
    println()
 
    return minimum(res), lcbtsvectoten(minimizer(res), m, n, d) 
end


"""
LC-BTS configuration sweep.

# Arguments
- `eta`: guess bts matrices with m levels and n pairs 
- `hamp`: Hamiltonian parameter 
- `swpiters`: number of sweep iterations
- `optiters`: number of optimization iterations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function swpconflcbts(eta::AbstractArray{T},
                      hamp::Float64,
                      swpiters::Int = 1000, 
                      optiters::Int = 5000, 
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    eps = 1e-08 
    en = randn(T)
    e0 = randn(T)
    etaten = zeros(d, n, m) 
    etaten .= copy(eta) 

    # Start the loop
    for j = 1:swpiters

        # Sweep
        for a = 1:d
            en, etamat = optonelcbts(etaten, a, hamp, optiters, hamtype, bound) 
            etaten[a, :, :] .= copy(etamat)
        end

        # Energy 
        ov = lcbtsoverlap(etaten) 
        en = lcbtshamoverlap(etaten, hamp, hamtype, bound) / ov 

        # Print
        println("--------")
        println("sweep: $j, en: $en")

        # Check convergence
        abs(en - e0) < eps && return en, etaten
        e0 = copy(en)

    end
 
    # Warning
    println("Warning: did not converge in $swpiters BTS sweep iterations")

    return en, etaten
end


"""
Optimize one BTS of the LC-BTS. 

# Arguments
- `eta`: guess bts matrices with m levels and n pairs 
- `btsind`: BTS index 
- `hamp`: Hamiltonian parameter 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function optonelcbts(eta::AbstractArray{T},
                     btsind::Int, 
                     hamp::Float64,
                     niters::Int = 1000, 
                     hamtype::Symbol=:RBCS,
                     bound::Symbol=:PBC) where T<:BlasFloat

    # Check
    1 <= btsind <= size(eta, 1) || error("Wrong BTS index in optonelcbts") 
    
    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    etavec = btpmattovec(eta[btsind, :, :])

    # Function
    function f(xvec)

        locten = zeros(d, n, m)
        locten[1:d-1, :, :] .= copy(eta[[1:btsind-1; btsind+1:d], :, :]) 
        locten[d, :, :] .= btpvectomat(xvec, m, n) 
        ov = lcbtsoverlap(locten) 
        en = lcbtshamoverlap(locten, hamp, hamtype, bound) 
        # println("elcbts with index $btsind:", ' ', en / ov)   

        return en / ov
    end
 
    # Optimize 
    @time res = optimize(f, etavec, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    return minimum(res), btpvectomat(minimizer(res), m, n) 
end


"""
<Psi|Psi> for a LC-BTS.

# Arguments
- `eta`: bts matrices of m levels and n pairs
- `cvec`: NOCI coefficients 
"""
function lcbtsoverlap(eta::AbstractArray{T},
                      cvec::AbstractVector{T}=ones(T, size(eta, 1))) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)

    # Scalar
    ov = zero(T)
    for i = 1:d
	ov += ( cvec[i]^2 * btsbraket(eta[i, :, :], eta[i, :, :]) )
    end	    
    for i = 1:d-1, j = i+1:d
	ov += ( 2 * conj(cvec[i]) * cvec[j] * btsbraket(eta[i, :, :], eta[j, :, :]) )
    end

    return ov 
end


"""
<Psi| H |Psi> for a LC-BTS.

# Arguments
- `eta`: bts matrices of m levels and n pairs
- `cvec`: NOCI coefficients 
- `hamp`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function lcbtshamoverlap(eta::AbstractArray{T},
                         hamp::Float64,
                         hamtype::Symbol=:RBCS,
                         bound::Symbol=:PBC,
                         cvec::AbstractVector{T}=ones(T, size(eta, 1))) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    
    # Scalar
    ov = zero(T)
    for i = 1:d
	ov += ( cvec[i]^2 * btsbrahamket(eta[i, :, :], eta[i, :, :], hamp, hamtype, bound) )
    end	    
    for i = 1:d-1, j = i+1:d
	ov += ( 2 * conj(cvec[i]) * cvec[j] * btsbrahamket(eta[i, :, :], eta[j, :, :], hamp, hamtype, bound) )
    end

    return ov 
end


"""
Pack BTS matrices into a vector.

# Argument
- `eta`: bts matrices of m levels and n pairs
"""
function lcbtstentovec(eta::AbstractArray{T}) where T<:BlasFloat

    # Initialize
    d = size(eta, 1)
    n = size(eta, 2)
    m = size(eta, 3)
    V = zeros(d * btpdim(m, n)) 

    # Trannsform
    a = 1
    b = 0
    for j = 1:d
        b = a + btpdim(m, n) - 1 
        V[a:b] = btpmattovec(eta[j, :, :]) 
        a = b+1
    end

    return V
end


"""
Unpack LC-BTS vector to BTS matrices. 

# Arguments
- `eta`: LC-BTS vector 
- `m`: Number of levels 
- `n`: Number of pairs 
- `d`: Number of BT states 
"""
function lcbtsvectoten(eta::AbstractVector{T},
                       m::Int,
                       n::Int,
                       d::Int) where T<:BlasFloat

    # Check 
    d*btpdim(m, n) == size(eta, 1) || error("Wrong dimension in lcbtsvectoten")  
    T3 = zeros(d, n, m) 

    # Transform
    a = 1
    b = 0
    for j = 1:d
        b = a + btpdim(m, n) - 1 
        T3[j, :, :] = btpvectomat(eta[a:b], m, n) 
        a = b+1
    end

    return T3 
end


