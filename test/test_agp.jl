using FermiPairing 
const fp = FermiPairing 
using Test

@testset "AGP for RBCS, m8n4, G = 0.30" begin

    # Initialize
    m = 8 
    n = 4 
    G = 0.30 

    # Tests
    etavec = pagpfrompcid(m, n, G)
    en1, eta = fp.pairingbfagp(m, n, G, etavec) 
    en2, eta = pairingagp(m, n, G, etavec) 
    @test isapprox(en1, 18.486123593452181, atol = 1e-8)
    @test isapprox(en1, en2, atol = 1e-8)

end


