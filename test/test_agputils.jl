using FermiPairing 
const fp = FermiPairing
using Test

@testset "AGP to HF" begin

    # Initialize
    m = rand(4:40)
    n = rand(2:m-1) 
    G = randn()
    eta = zeros(Float64, m) 

    # Tests
    eta[1:n] = ones(Float64, n)
    m1 = agpoverlap(eta, eta, n) 
    h1 = agphamoverlap(eta, eta, n, G, :RBCS) 
    h2 = agphamoverlap(eta, eta, n, G, :XXZ, :OBC) 
    h3 = agphamoverlap(eta, eta, n, G, :XXZ, :PBC) 
    en1 = h1 / m1
    en2 = h2 / m1
    en3 = h3 / m1
    en4 = pairinghf(m, n, G, :RBCS)
    en5 = pairinghf(m, n, G, :XXZ, :OBC)
    en6 = pairinghf(m, n, G, :XXZ, :PBC)
    @test isapprox(en1, en4, atol = 1e-12)
    @test isapprox(en2, en5, atol = 1e-12)
    @test isapprox(en3, en6, atol = 1e-12)

end


@testset "AGP overlaps" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1) 
    G = randn()
    eta1 = randn(m)
    eta2 = randn(m)

    # Tests
    m1 = agpoverlap(eta1, eta2, n) 
    h1 = agphamoverlap(eta1, eta2, n, G, :RBCS) 
    h2 = agphamoverlap(eta1, eta2, n, G, :XXZ, :OBC) 
    h3 = agphamoverlap(eta1, eta2, n, G, :XXZ, :PBC) 
    h4, m2 = fp.bfagpoverlaps(eta1, eta2, n, G, :RBCS) 
    h5, m2 = fp.bfagpoverlaps(eta1, eta2, n, G, :XXZ, :OBC) 
    h6, m2 = fp.bfagpoverlaps(eta1, eta2, n, G, :XXZ, :PBC) 
    @test isapprox(m1, m2, atol = 1e-12)
    @test isapprox(h1, h4, atol = 1e-12)
    @test isapprox(h2, h5, atol = 1e-12)
    @test isapprox(h3, h6, atol = 1e-12)

end

