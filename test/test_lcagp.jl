using FermiPairing 
const fp = FermiPairing 
using Test

@testset "LCAGP and AGP" begin

    # Initialize
    m = rand(4:40)
    n = rand(2:m-1) 
    G = randn()
    eta1 = normalizeagp!(randn(m), n)
    eta2 = normalizeagp!(randn(m), n)

    # Tests
    en1 = agphamoverlap(eta1, eta1, n, G, :RBCS) / agpoverlap(eta1, eta1, n)
    en2 = agphamoverlap(eta1, eta1, n, G, :XXZ, :OBC) / agpoverlap(eta1, eta1, n)
    en3 = agphamoverlap(eta1, eta1, n, G, :XXZ, :PBC) / agpoverlap(eta1, eta1, n)
    en4, evec, zm = pairinglcagp_gen(transpose(hcat(eta1, eta2)), n, G, :RBCS) 
    en5, evec, zm = pairinglcagp_gen(transpose(hcat(eta1, eta2)), n, G, :XXZ, :OBC) 
    en6, evec, zm = pairinglcagp_gen(transpose(hcat(eta1, eta2)), n, G, :XXZ, :PBC) 
    @test en4 < en1 
    @test en5 < en2 
    @test en6 < en3 

end


@testset "J1-CI" begin

    # Initialize
    m = rand(4:10)
    n = 1 
    G = randn() 
    eta = normalizeagp!(randn(m), n)
    pivot = 0.0

    # Tests
    e1, evec, zm = pairinglcagp_pivot(eta, n, pivot, 1, G, :RBCS) 
    e2, V = pairingfci(m, n, G, :RBCS)
    e3, evec, zm = pairinglcagp_pivot(eta, n, pivot, 1, G, :XXZ, :OBC) 
    e4, V = pairingfci(m, n, G, :XXZ, :OBC)
    e5, evec, zm = pairinglcagp_pivot(eta, n, pivot, 1, G, :XXZ, :PBC) 
    e6, V = pairingfci(m, n, G, :XXZ, :PBC)
    @test isapprox(e1, e2, atol = 1e-8)
    @test isapprox(e3, e4, atol = 1e-8)
    @test isapprox(e5, e6, atol = 1e-8)

end


@testset "J2-CI" begin

    # Initialize
    m = rand(4:10)
    n = 2 
    G = randn() 
    eta = normalizeagp!(randn(m), n)
    pivot = 0.0

    # Tests
    e1, evec, zm = pairinglcagp_pivot(eta, n, pivot, 2, G, :RBCS) 
    e2, V = pairingfci(m, n, G, :RBCS)
    e3, evec, zm = pairinglcagp_pivot(eta, n, pivot, 2, G, :XXZ, :OBC) 
    e4, V = pairingfci(m, n, G, :XXZ, :OBC)
    e5, evec, zm = pairinglcagp_pivot(eta, n, pivot, 2, G, :XXZ, :PBC) 
    e6, V = pairingfci(m, n, G, :XXZ, :PBC)
    @test isapprox(e1, e2, atol = 1e-8)
    @test isapprox(e3, e4, atol = 1e-8)
    @test isapprox(e5, e6, atol = 1e-8)

end
