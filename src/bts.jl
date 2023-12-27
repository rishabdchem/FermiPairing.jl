using LinearAlgebra 
using Optim
using Optim: converged, iterations, minimum, minimizer

import LinearAlgebra.BLAS: BlasFloat
export btsmin, btsfromagp 

#==================================
Solve single BTS.
==================================#

"""
ebts, etamat <- optimize trial BTS ground state energy.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamp`: Hamiltonian parameter 
- `eta`: guess bts 
- `niters`: number of itrations
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: boundary condition with PBC as default 
- `eps`: pairing levels 
"""
function btsmin(m::Int,
                n::Int, 
                hamp::Float64,
                eta::AbstractVector{T},
                niters::Int = 1000, 
                hamtype::Symbol=:RBCS,
                bound::Symbol=:PBC, 
                eps::Vector{Int}=Vector(1:size(eta, 1))) where T<:BlasFloat

    # Check
    size(eta, 1) == btpdim(m, n) || error("Wrong BTS vector dimension bfbtsmin")

    # Function
    function f(xvec)

        ov = btsoverlap(btpvectomat(xvec, m, n)) 
        en = btshamoverlap(btpvectomat(xvec, m, n), hamp, hamtype, bound, eps) 
        # println("ebts inside:", ' ', en / ov)   

        return en / ov
    end
 
    # Print
    println("----------------")
    println("Initial BTS:")
    display(btpvectomat(eta, m, n))
    println()
    
    # Optimize 
    @time res = optimize(f, eta, ConjugateGradient(), 
                         Optim.Options(f_tol = 1e-12,
                                       iterations = niters))

    # Check convergence
    converged(res) || println("Warning: did not converge after $(iterations(res)) iterations")

    # Print
    println("------------BTS details-------------")
    display(res)
    println()

    return minimum(res), btpvectomat(minimizer(res), m, n) 
end


"""
Guess BTS from AGP. 

# Arguments
- `eta`: AGP of m levels 
- `n`: number of pairs
- `sympar`: symmetry breaking paramter  
"""
function btsfromagp(eta::AbstractVector{T},
                    n::Int,
                    sympar::Float64=1e-06) where T<:BlasFloat

    # Initialize
    m = size(eta, 1)
    G = zeros(T, n, m)
    G[1, 1] = copy(eta[1])

    # Guess 
    for p = 2:m, j = btprowmin(m, n, p):btprowmax(n, p) 
        G[j, p] = copy(eta[p]) 
        j == btprowmin(m, n, p) && continue
        G[j, p] += (sympar * randn(T)) 
    end

    return btpmattovec(G) 
end


"""
Guess BT states from HF-CID.

# Arguments
- `m`: number of levels 
- `n`: number of pairs
- `hamparam`: Hamiltonian parameter 
- `hamtype`: Hamiltonian type with RBCS as default 
- `bound`: Boundary condition with PBC as default 
- `nstates`: number of states 
"""
function btstatesfromcid(m::Int,
                         n::Int,
                         hamparam::Float64,
                         hamtype::Symbol=:RBCS,
                         bound::Symbol=:PBC,
                         nstates::Int=1) 

    # Check 
    n > m && error("n > m in btstatesfromcid") 
    1 <= nstate <= min(n, m - n) || error("Wrong number of BTS in btstatesfromcid") 

    # Get AGPs
    etamat = pagpfrompcid(m, n, hamparam, hamtype, bound, nstates)   

    # Get BT states
    btsten = zeros(Float64, nstate, n, m)
    for j = 1:nstate
        btsten[j, :, :] = btpvectomat(btsfromagp(etamat[j, :]), m, n)
    end 

    return btsten 
end


