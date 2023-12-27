using FermiPairing 
const fp = FermiPairing 
using Test

@testset "Res-AGP" begin

    # Initialize
    m = 4 
    n = 2 
    G = randn()
    eta = zeros(Float64, 2, m)
    eta[1, :] = normalizeagp!(randn(m), n)
    eta[2, :] = normalizeagp!(randn(m), n)

    # Tests
    en1 = agphamoverlap(eta[1, :], eta[1, :], n, G, :RBCS) / agpoverlap(eta[1, :], eta[1, :], n)
    en2 = agphamoverlap(eta[1, :], eta[1, :], n, G, :XXZ, :OBC) / agpoverlap(eta[1, :], eta[1, :], n)
    en3 = agphamoverlap(eta[1, :], eta[1, :], n, G, :XXZ, :PBC) / agpoverlap(eta[1, :], eta[1, :], n)
    en4, etamat = pairingresagp(eta, n, G, 1000, :RBCS) 
    en5, etamat = pairingresagp(eta, n, G, 1000, :XXZ, :OBC) 
    en6, etamat = pairingresagp(eta, n, G, 1000, :XXZ, :PBC) 
    @test en4 < en1 
    @test en5 < en2 
    @test en6 < en3 

end


@testset "Sweep-LC-AGP" begin

    # Initialize
    m = 4 
    n = 2
    G = randn()
    eta = zeros(Float64, 2, m)
    eta[1, :] = normalizeagp!(randn(m), n)
    eta[2, :] = normalizeagp!(randn(m), n)

    # Tests
    en1 = agphamoverlap(eta[1, :], eta[1, :], n, G, :RBCS) / agpoverlap(eta[1, :], eta[1, :], n)
    en2 = agphamoverlap(eta[1, :], eta[1, :], n, G, :XXZ, :OBC) / agpoverlap(eta[1, :], eta[1, :], n)
    en3 = agphamoverlap(eta[1, :], eta[1, :], n, G, :XXZ, :PBC) / agpoverlap(eta[1, :], eta[1, :], n)
    en4, etamat = pairingswpconflcagp(eta, n, G, 10, 1000, :RBCS) 
    en5, etamat = pairingswpconflcagp(eta, n, G, 10, 1000, :XXZ, :OBC) 
    en6, etamat = pairingswpconflcagp(eta, n, G, 10, 1000, :XXZ, :PBC) 
    @test en4 < en1 
    @test en5 < en2 
    @test en6 < en3 

end


