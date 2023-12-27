using FermiPairing 
const fp = FermiPairing 
using Test

@testset "APIG to AGP" begin

    # Initialize
    m = 4 
    n = 2 
    G = randn()
    eta = randn(m)
    etamat = zeros(Float64, n, m) 

    # Tests
    for i = 1:n
        etamat[i, :] .= copy(eta)
    end
    m1 = agpoverlap(eta, eta, n) 
    h1 = agphamoverlap(eta, eta, n, G, :RBCS) 
    h2, m2 = fp.apigoverlaps(etamat, G, :RBCS) 
    h3 = agphamoverlap(eta, eta, n, G, :XXZ, :PBC) 
    h4, m2 = fp.apigoverlaps(etamat, G, :XXZ, :PBC) 
    h5 = agphamoverlap(eta, eta, n, G, :XXZ, :OBC) 
    h6, m2 = fp.apigoverlaps(etamat, G, :XXZ, :OBC) 
    @test isapprox(h1/m1, h2/m2, atol = 1e-08)
    @test isapprox(h3/m1, h4/m2, atol = 1e-08)
    @test isapprox(h5/m1, h6/m2, atol = 1e-08)

end


@testset "Compare APIG and AGP" begin

    # Initialize
    m = 5 
    n = 3 
    G = randn()
    niters = 1000
    println("m: $m, n: $n, GorD: $G")
 
    # Tests
    etavec = pagpfrompcid(m, n, G, :RBCS)
    en1, eta = pairingagp(m, n, G, etavec, niters, :RBCS) 
    etaguess = papigfrompagp(eta, n) 
    en2, etamat = pairingapig(m, n, G, etaguess, niters, :RBCS) 
    etavec = pagpfrompcid(m, n, G, :XXZ, :PBC)
    en3, eta = pairingagp(m, n, G, etavec, niters, :XXZ, :PBC) 
    etaguess = papigfrompagp(eta, n) 
    en4, etamat = pairingapig(m, n, G, etaguess, niters, :XXZ, :PBC) 
    etavec = pagpfrompcid(m, n, G, :XXZ, :OBC)
    en5, eta = pairingagp(m, n, G, etavec, niters, :XXZ, :OBC) 
    etaguess = papigfrompagp(eta, n) 
    en6, etamat = pairingapig(m, n, G, etaguess, niters, :XXZ, :OBC) 
    @test en2 < en1 
    @test en4 < en3 
    @test en6 < en5 

end


