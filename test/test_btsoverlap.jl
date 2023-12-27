using FermiPairing 
const fp = FermiPairing
using Test

@testset "BTS normalization" begin

    # Initialize
    m = rand(4:40)
    n = rand(2:m-1)
    G = randn() 
    B = randn(n, m)

    # Tests
    m1 = btsoverlap(B)
    h1 = btshamoverlap(B, G, :RBCS)
    h2 = btshamoverlap(B, G, :XXZ, :PBC)
    h3 = btshamoverlap(B, G, :XXZ, :OBC)
    m2 = btsoverlap(normalizebts!(B))
    h4 = btshamoverlap(normalizebts!(B), G, :RBCS)
    h5 = btshamoverlap(normalizebts!(B), G, :XXZ, :PBC)
    h6 = btshamoverlap(normalizebts!(B), G, :XXZ, :OBC)
    en1 = h1/m1
    en2 = h2/m1
    en3 = h3/m1
    en4 = h4/m2
    en5 = h5/m2
    en6 = h6/m2
    println("m: $m, n: $n")
    println("New overlap: $m2")
    @test isapprox(1.0, m2, atol=1e-12) 
    @test isapprox(en1, en4, atol=1e-08) 
    @test isapprox(en2, en5, atol=1e-08) 
    @test isapprox(en3, en6, atol=1e-08) 

end


@testset "BTS overlaps" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    G = randn()
    B = randn(n, m)

    # Tests
    m1 = btsoverlap(B)   
    h1 = btshamoverlap(B, G)   
    h2 = btshamoverlap(B, G, :XXZ, :PBC)   
    h3 = btshamoverlap(B, G, :XXZ, :OBC)   
    h4, m2 = fp.bfbtsoverlaps(B, G)   
    h5, m2 = fp.bfbtsoverlaps(B, G, :XXZ, :PBC)   
    h6, m2 = fp.bfbtsoverlaps(B, G, :XXZ, :OBC)   
    @test isapprox(m1, m2, atol=1e-08) 
    @test isapprox(h1, h4, atol=1e-08) 
    @test isapprox(h2, h5, atol=1e-08) 
    @test isapprox(h3, h6, atol=1e-08) 

end


@testset "BTS transition overlaps" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    G = randn()
    B = randn(n, m)
    C = randn(n, m)

    # Tests
    m1 = btsbraket(B, C)   
    h1 = btsbrahamket(B, C, G)   
    h2 = btsbrahamket(B, C, G, :XXZ, :PBC)   
    h3 = btsbrahamket(B, C, G, :XXZ, :OBC)   
    h4, m2 = fp.bfbtstransitions(B, C, G)   
    h5, m2 = fp.bfbtstransitions(B, C, G, :XXZ, :PBC)   
    h6, m2 = fp.bfbtstransitions(B, C, G, :XXZ, :OBC)   
    @test isapprox(m1, m2, atol=1e-08) 
    @test isapprox(h1, h4, atol=1e-08) 
    @test isapprox(h2, h5, atol=1e-08) 
    @test isapprox(h3, h6, atol=1e-08) 

end
