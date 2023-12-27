using FermiPairing 
const fp = FermiPairing
using Combinatorics 
using Test

@testset "BTS TDMs" begin

    # Initialize
    m = rand(7:10)
    n = rand(4:m-1)
    d = binomial(m, n)
    B = randn(n, m)
    C = randn(n, m)
    U = ones(Float64, d)
    V = ones(Float64, d)
    p = rand(2:m-1)
    q = rand([2:p-1; p+1:m-1]) 
    println("p is $p and q is $q")

    # Get combinations
    X = Vector(1:m)
    comb = combinations(X, n)

    # Coefficients 
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            U[mu] *= B[j, pind[j]]
            V[mu] *= C[j, pind[j]]
        end    
    end 

    # Tests
    m1 = fp.btstz11(B, C, p)   
    m2 = fp.btstz22(B, C, (p, q,))   
    m3 = fp.btstz22(B, C, (p, 1,))   
    m4 = fp.btstz22(B, C, (p, n,))   
    m5 = fp.btstz22(B, C, (p, m,))   
    m6 = fp.btstz22(B, C, (q, 1,))   
    m7 = fp.btstz22(B, C, (q, n,))   
    m8 = fp.btstz22(B, C, (q, m,))  
    s1 = fp.docibranumket(U, V, m, n, p) 
    s2 = fp.docibranumnumket(U, V, m, n, (p, q,)) 
    s3 = fp.docibranumnumket(U, V, m, n, (p, 1,)) 
    s4 = fp.docibranumnumket(U, V, m, n, (p, n,)) 
    s5 = fp.docibranumnumket(U, V, m, n, (p, m,)) 
    s6 = fp.docibranumnumket(U, V, m, n, (q, 1,)) 
    s7 = fp.docibranumnumket(U, V, m, n, (q, n,)) 
    s8 = fp.docibranumnumket(U, V, m, n, (q, m,)) 
    @test isapprox(m1, s1, atol=1e-12) 
    @test isapprox(m2, s2, atol=1e-12) 
    @test isapprox(m3, s3, atol=1e-12) 
    @test isapprox(m4, s4, atol=1e-12) 
    @test isapprox(m5, s5, atol=1e-12) 
    @test isapprox(m6, s6, atol=1e-12) 
    @test isapprox(m7, s7, atol=1e-12) 
    @test isapprox(m8, s8, atol=1e-12) 

    # Tests
    t1 = fp.btstz02(B, C, (p, q,))   
    t2 = fp.btstz02(B, C, (q, p,))   
    t3 = fp.btstz02(B, C, (p, 1,))   
    t4 = fp.btstz02(B, C, (1, p,))   
    t5 = fp.btstz02(B, C, (p, n,))   
    t6 = fp.btstz02(B, C, (n, p,))   
    t7 = fp.btstz02(B, C, (p, m,))   
    t8 = fp.btstz02(B, C, (m, p,))   
    t9 = fp.btstz02(B, C, (q, 1,))   
    t10 = fp.btstz02(B, C, (1, q,))   
    t11 = fp.btstz02(B, C, (q, n,))   
    t12 = fp.btstz02(B, C, (n, q,))   
    t13 = fp.btstz02(B, C, (q, m,))   
    t14 = fp.btstz02(B, C, (m, q,))   
    u1 = fp.docibrapdagpket(U, V, m, n, (p, q,)) 
    u2 = fp.docibrapdagpket(U, V, m, n, (q, p,)) 
    u3 = fp.docibrapdagpket(U, V, m, n, (p, 1,)) 
    u4 = fp.docibrapdagpket(U, V, m, n, (1, p,)) 
    u5 = fp.docibrapdagpket(U, V, m, n, (p, n,)) 
    u6 = fp.docibrapdagpket(U, V, m, n, (n, p,)) 
    u7 = fp.docibrapdagpket(U, V, m, n, (p, m,)) 
    u8 = fp.docibrapdagpket(U, V, m, n, (m, p,)) 
    u9 = fp.docibrapdagpket(U, V, m, n, (q, 1,)) 
    u10 = fp.docibrapdagpket(U, V, m, n, (1, q,)) 
    u11 = fp.docibrapdagpket(U, V, m, n, (q, n,)) 
    u12 = fp.docibrapdagpket(U, V, m, n, (n, q,)) 
    u13 = fp.docibrapdagpket(U, V, m, n, (q, m,)) 
    u14 = fp.docibrapdagpket(U, V, m, n, (m, q,)) 
    @test isapprox(t1, u1, atol=1e-12) 
    @test isapprox(t2, u2, atol=1e-12) 
    @test isapprox(t3, u3, atol=1e-12) 
    @test isapprox(t4, u4, atol=1e-12) 
    @test isapprox(t5, u5, atol=1e-12) 
    @test isapprox(t6, u6, atol=1e-12) 
    @test isapprox(t7, u7, atol=1e-12) 
    @test isapprox(t8, u8, atol=1e-12) 
    @test isapprox(t9, u9, atol=1e-12) 
    @test isapprox(t10, u10, atol=1e-12) 
    @test isapprox(t11, u11, atol=1e-12) 
    @test isapprox(t12, u12, atol=1e-12) 
    @test isapprox(t13, u13, atol=1e-12) 
    @test isapprox(t14, u14, atol=1e-12) 

end


@testset "Number operator on BTS" begin

    # Initialize
    m = rand(7:10)
    n = rand(4:m-1)
    d = binomial(m, n)
    B = randn(n, m)
    C = copy(B) 
    U = ones(Float64, d)
    V = ones(Float64, d)

    # Number operator 
    p = rand(1:m)
    beta = 2 
    C[:, p] .*= (1 + (2 / beta)) 

    # Coefficients 
    comb = combinations(Vector(1:m), n)
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            U[mu] *= B[j, pind[j]]
            V[mu] *= C[j, pind[j]]
        end    
    end 
    V .*= beta

    # Test
    W = 2 .* actopdoci(U, m, n, (p,), (p,)) 
    W .+= (beta .* U) 
    @test isapprox(V, W, atol=1e-12) 
    
end


@testset "BTS number operator overlaps" begin

    # Initialize
    m = rand(7:10)
    n = rand(4:m-1)
    d = binomial(m, n)
    B = randn(n, m)
    C = randn(n, m)
    U = ones(Float64, d)
    V = ones(Float64, d)

    # Different indices assumed 
    p = rand(1:m-3)
    q = p+1
    r = q+1
    s = r+1
    println("m: $m, n: $n")
    println("p, q, r, s: $p, $q, $r, $s")

    # Coefficients 
    comb = combinations(Vector(1:m), n)
    mu = 0
    for pind in comb
        mu += 1
        for j = 1:n    
            U[mu] *= B[j, pind[j]]
            V[mu] *= C[j, pind[j]]
        end    
    end 

    # Reference values
    V1 = 2 .* actopdoci(U, m, n, (p,), (p,)) 
    V2 = 2 .* actopdoci(V, m, n, (s,), (s,)) 
    r1 = fp.docibrapdagpket(V1, V2, m, n, (q, r,)) 
    r2 = fp.docibrapdagpket(V1, V, m, n, (q, r,)) 
    r3 = fp.docibrapdagpket(U, V2, m, n, (q, r,)) 
    r4 = fp.docibranumnumket(V1, V2, m, n, (q, r,)) 

    # Tests
    t1 = fp.btsnumz02num(B, C, (p, q, r, s,))
    t2 = fp.btsnumz02(B, C, (p, q, r,))
    t3 = fp.btsz02num(B, C, (q, r, s,))
    t4 = fp.btstz44(B, C, (p, q, r, s,))
    @test isapprox(r1, t1, atol=1e-12) 
    @test isapprox(r2, t2, atol=1e-12) 
    @test isapprox(r3, t3, atol=1e-12) 
    @test isapprox(r4, t4, atol=1e-12) 
    
end

