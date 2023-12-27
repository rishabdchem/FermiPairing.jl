using FermiPairing 
using Test

@testset "LC-BTS and BTS" begin

    # Initialize
    m = rand(4:40)
    n = rand(2:m-1) 
    G = randn()
    eta1 = normalizebts!(randn(n, m))
    eta2 = normalizebts!(randn(n, m))
    etaten = zeros(2, n, m) 
    etaten[1, :, :] .= copy(eta1)
    etaten[2, :, :] .= copy(eta2)

    # Tests
    en1 = btshamoverlap(eta1, G, :RBCS) / btsoverlap(eta1)
    en2 = btshamoverlap(eta1, G, :XXZ, :OBC) / btsoverlap(eta1)
    en3 = btshamoverlap(eta1, G, :XXZ, :PBC) / btsoverlap(eta1)
    en4, evec, zm = genlcbts(etaten, G, :RBCS) 
    en5, evec, zm = genlcbts(etaten, G, :XXZ, :OBC) 
    en6, evec, zm = genlcbts(etaten, G, :XXZ, :PBC) 
    @test en4 < en1 
    @test en5 < en2 
    @test en6 < en3 

end

