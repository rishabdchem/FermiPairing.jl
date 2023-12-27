using FermiPairing 
using Test

@testset "Compare BTS guesses" begin

    # Initialize
    m = rand(4:10) 
    n = rand(2:m-1)
    G = randn()
    println("m: $m, n: $n")

    # Tests
    etavec = pagpfrompcid(m, n, G, :RBCS) 
    eagp, eta = pairingagp(m, n, G, etavec, 5000)
    btsvec = btsfromagp(eta, n)
    e1, btsmat = btsmin(m, n, G, btsvec, 5000)
    btsvec = btsfromagp(etavec, n)
    e2, btsmat = btsmin(m, n, G, btsvec, 5000)
    @test isapprox(e1, e2, atol = 1e-8)

end


@testset "BTS for RBCS, m10n5" begin

    # Initialize
    m = 10 
    n = 5

    # Test
    G = 0.30 
    eta = pagpfrompcid(m, n, G, :RBCS) 
    btsvec = btsfromagp(eta, n)
    en, btsmat = btsmin(m, n, G, btsvec, 5000)
    @test isapprox(en, 28.072303343960094, atol = 1e-6)

    # Test
    G = -0.50 
    eta = pagpfrompcid(m, n, G, :RBCS) 
    btsvec = btsfromagp(eta, n)
    en, btsmat = btsmin(m, n, G, btsvec, 5000)
    @test isapprox(en, 31.985730172933167, atol = 1e-6)

end


