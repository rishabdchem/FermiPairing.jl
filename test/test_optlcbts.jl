using FermiPairing 
const fp = FermiPairing
using Test

@testset "Res-BTS and BTS" begin

    # Initialize
    m = 8 
    n = 4 
    G = randn()

    # BTS 
    etavec = pagpfrompcid(m, n, G, :RBCS) 
    eagp, eta = pairingagp(m, n, G, etavec, 5000)
    btsvec = btsfromagp(eta, n)
    en1, btsmat = btsmin(m, n, G, btsvec, 5000)

    # Res-BTS
    eta1 = normalizebts!(btsmat)
    eta2 = normalizebts!(randn(n, m))
    etaten = zeros(2, n, m) 
    etaten[1, :, :] = copy(eta1)
    etaten[2, :, :] = copy(eta2)

    # Tests
    en2, etanew = resbts(etaten, G, 1000, :RBCS) 
    @test en2 < en1 

end

