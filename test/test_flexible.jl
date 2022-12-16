using FermiPairing 
const fp = FermiPairing 
using Test

@testset "HF for RBCS" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    G = randn()
    V = zeros(Float64, binomial(m, n))

    # Tests
    V[1] = 1.0
    en = 0.0
    for i = 1:n
        en += 2i - G
    end  
    @test isapprox(1.0, docioverlap(V, V, m, n), atol=1e-12) 
    @test isapprox(en, docihamoverlap(V, V, m, n, G), atol=1e-12) 

end


@testset "Flexible RDMs" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    d = binomial(m, n)
    U = randn(d) 
    V = randn(d) 
    G = randn()

    # Tests
    h1 = 0.0
    for p = 1:m
        h1 += (p - G/2) * fp.docibranumket(U, V, m, n, p) 
    end 
    for p = 1:m, q = 1:m
        p == q && continue
        h1 -= G * fp.docibrapdagpket(U, V, m, n, (p, q,)) 
    end
    h2 = docihamoverlap(U, V, m, n, G)
    @test isapprox(h1, h2, atol=1e-12) 

end


@testset "Flexible operator" begin

    # Initialize
    m = rand(4:10)
    n = rand(2:m-1)
    d = binomial(m, n)
    U = randn(d) 
    V = randn(d) 
    X = randn(d) 
    G = randn()

    # Tests
    h1 = 0.0
    for p = 1:m
        X = actopdoci(V, m, n, (p,), (p,))
        h1 += (2p - G) * docioverlap(U, X, m, n) 
    end 
    for p = 1:m, q = 1:m
        p == q && continue
        X = actopdoci(V, m, n, (p,), (q,))
        h1 -= G * docioverlap(U, X, m, n) 
    end
    h2 = docihamoverlap(U, V, m, n, G)
    @test isapprox(h1, h2, atol=1e-12) 

end


