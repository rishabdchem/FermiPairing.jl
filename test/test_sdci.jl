using FermiPairing 
const fp = FermiPairing 
using Test

@testset "HF" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    G = randn()
    V = zeros(Float64, binomial(m, n))

    # Tests
    V[1] = 1.0
    e1 = pairinghf(m, n, G, :RBCS)
    e2 = pairinghf(m, n, G, :XXZ, :OBC)
    e3 = pairinghf(m, n, G, :XXZ, :PBC)
    @test isapprox(e1, docihamoverlap(V, V, m, n, G), atol=1e-12) 
    @test isapprox(e2, docihamoverlap(V, V, m, n, G, :XXZ, :OBC), atol=1e-12) 
    @test isapprox(e3, docihamoverlap(V, V, m, n, G, :XXZ, :PBC), atol=1e-12) 

end


@testset "CID" begin

    # Initialize
    m = rand(4:16)
    n = 1 
    G = randn()

    # Tests
    e1, V1 = pairingfci(m, n, G)
    e2, V1 = pairingcid(m, n, G)
    e3, V2 = pairingfci(m, n, G, :XXZ, :OBC)
    e4, V2 = pairingcid(m, n, G, :XXZ, :OBC)
    e5, V3 = pairingfci(m, n, G, :XXZ, :PBC)
    e6, V3 = pairingcid(m, n, G, :XXZ, :PBC)
    @test isapprox(e1, e2, atol=1e-8) 
    @test isapprox(e3, e4, atol=1e-8) 
    @test isapprox(e5, e6, atol=1e-8) 

end


@testset "FCI" begin

    # Initialize
    m = rand(4:16)
    n = rand(2:m-1)
    G = randn()

    # Tests
    e1, V1 = pairingfci(m, n, G)
    e2, V2 = pairingfci(m, n, G, :XXZ, :OBC)
    e3, V3 = pairingfci(m, n, G, :XXZ, :PBC)
    @test isapprox(e1, docihamoverlap(V1, V1, m, n, G), atol=1e-12) 
    @test isapprox(e2, docihamoverlap(V2, V2, m, n, G, :XXZ, :OBC), atol=1e-12) 
    @test isapprox(e3, docihamoverlap(V3, V3, m, n, G, :XXZ, :PBC), atol=1e-12) 

end
