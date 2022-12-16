import LinearAlgebra.BLAS: BlasFloat
using Optim 
using Optim: converged, iterations, minimum, minimizer
using Optim: only_fg!, only_fgh!

export pairingagp, pagpfrompcid

#==================================
Single AGP solver. 
==================================#

"""
Optimize trial AGP ground state energy.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `eta`: guess AGP 
- `niters`: number of iterations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairingagp(m::Int,
                    n::Int, 
                    hamp::Float64,
                    eta::AbstractVector{T},
                    niters::Int = 1000, 
                    hamtype::Symbol=:RBCS,
                    bound::Symbol=:PBC) where T<:BlasFloat

    # Function
    function f(xvec)

        ov = agpoverlap(xvec, xvec, n)
        en = agphamoverlap(xvec, xvec, n, hamp, hamtype, bound) 
        # println("eagp inside:", ' ', en / ov)   
 
        return en / ov 
    end

    # Print
    println("----------------")
    println("Initial AGP:")
    display(eta)
    println()
    
    # Optimize 
    @time res = optimize(f, eta, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || error("Did not converge after $(iterations(res)) iterations")

    # Print
    println("------------AGP details-------------")
    display(res)
    println()

    return minimum(res), minimizer(res) 
end


"""
Guess AGP from particle-hole doubles amplitude.
See DOI: 10.1103/PhysRevB.93.125124.

We get T2 from CID and after doing SVD

  T(a, i) = sum(mu) S(mu) U(a, mu) V(i, mu), 

take the modes corresponding to the largest singular value.

Returns real-valued geminal coefficients.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: Boundary condition with PBC as default 
- `agpnum': Number of AGPs to get  
- `enret': Return CID energy or not 
"""
function pagpfrompcid(m::Int,
                      n::Int, 
                      hamparam::Float64,
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC,
                      agpnum::Int=1,  
                      enret::Bool=false)

    # Check 
    n > m && error("n > m in pagpfrompcid") 
    agpnum > min(n, m - n) && error("Asking for too many AGPs in pagpfrompcid") 
    
    # CID
    en, C = pairingcid(m, n, hamparam, hamtype, bound) 
    n == 1 && return C
    n == m - 1 && return C

    # Initialize
    T = zeros(Float64, m-n, n)
    if agpnum == 1
        V = zeros(Float64, m)
    else
        V = zeros(Float64, agpnum, m)
    end 

    # Unpack C
    ai = 1
    for i = 1:n, a = 1:m-n 
        ai += 1 
        T[a, i] = C[ai]
    end

    # SVD
    F = svd(T)
    
    # Guess 
    if agpnum == 1
        V[1:n] = 1.0 ./ F.V[:, 1]
        V[n+1:m] = F.U[:, 1]
    else   
        for j = 1:agpnum
            V[j, 1:n] = 1.0 ./ F.V[:, j]
            V[j, n+1:m] = F.U[:, j]
        end
    end

    enret == true && return en, V
    return V 
end


"""
eagp, eta <- optimize trial AGP ground state energy.
Brute force computation of energies.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `eta`: guess AGP 
- `niters`: number of iterations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
"""
function pairingbfagp(m::Int,
                      n::Int, 
                      hamp::Float64,
                      eta::AbstractVector{T},
                      niters::Int = 1000, 
                      hamtype::Symbol=:RBCS,
                      bound::Symbol=:PBC) where T<:BlasFloat

    # Function
    function f(xvec)

        en, ov = bfagpoverlaps(xvec, xvec, n, hamp, hamtype, bound) 
        # println("eagp inside:", ' ', en / ov)   
 
        return en / ov 
    end

    # Print
    println("----------------")
    println("Initial AGP:")
    display(eta)
    println()
    
    # Optimize 
    @time res = optimize(f, eta, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || error("Did not converge after $(iterations(res)) iterations")

    # Print
    println("------------AGP details-------------")
    display(res)
    println()

    return minimum(res), minimizer(res) 
end

