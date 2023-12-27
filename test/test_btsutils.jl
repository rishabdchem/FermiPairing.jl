using FermiPairing 
const fp = FermiPairing
using Combinatorics 
using Test

@testset "BTS RDMs" begin

    # Initialize
    m = rand(7:10)
    n = rand(4:m-1)
    d = binomial(m, n)
    B = randn(n, m)
    U = ones(Float64, d)
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
        end    
    end 

    # Tests
    m1 = fp.btsz11(B, p)   
    m2 = fp.btsz22(B, (p, q,))   
    m3 = fp.btsz22(B, (p, 1,))   
    m4 = fp.btsz22(B, (p, n,))   
    m5 = fp.btsz22(B, (p, m,))   
    m6 = fp.btsz22(B, (q, 1,))   
    m7 = fp.btsz22(B, (q, n,))   
    m8 = fp.btsz22(B, (q, m,))  
    s1 = fp.docibranumket(U, U, m, n, p) 
    s2 = fp.docibranumnumket(U, U, m, n, (p, q,)) 
    s3 = fp.docibranumnumket(U, U, m, n, (p, 1,)) 
    s4 = fp.docibranumnumket(U, U, m, n, (p, n,)) 
    s5 = fp.docibranumnumket(U, U, m, n, (p, m,)) 
    s6 = fp.docibranumnumket(U, U, m, n, (q, 1,)) 
    s7 = fp.docibranumnumket(U, U, m, n, (q, n,)) 
    s8 = fp.docibranumnumket(U, U, m, n, (q, m,)) 
    @test isapprox(m1, s1, atol=1e-12) 
    @test isapprox(m2, s2, atol=1e-12) 
    @test isapprox(m3, s3, atol=1e-12) 
    @test isapprox(m4, s4, atol=1e-12) 
    @test isapprox(m5, s5, atol=1e-12) 
    @test isapprox(m6, s6, atol=1e-12) 
    @test isapprox(m7, s7, atol=1e-12) 
    @test isapprox(m8, s8, atol=1e-12) 

    # Tests
    t1 = fp.btsz02(B, (p, q,))   
    t2 = fp.btsz02(B, (q, p,))   
    t3 = fp.btsz02(B, (p, 1,))   
    t4 = fp.btsz02(B, (1, p,))   
    t5 = fp.btsz02(B, (p, n,))   
    t6 = fp.btsz02(B, (n, p,))   
    t7 = fp.btsz02(B, (p, m,))   
    t8 = fp.btsz02(B, (m, p,))   
    t9 = fp.btsz02(B, (q, 1,))   
    t10 = fp.btsz02(B, (1, q,))   
    t11 = fp.btsz02(B, (q, n,))   
    t12 = fp.btsz02(B, (n, q,))   
    t13 = fp.btsz02(B, (q, m,))   
    t14 = fp.btsz02(B, (m, q,))   
    u1 = fp.docibrapdagpket(U, U, m, n, (p, q,)) 
    u2 = fp.docibrapdagpket(U, U, m, n, (q, p,)) 
    u3 = fp.docibrapdagpket(U, U, m, n, (p, 1,)) 
    u4 = fp.docibrapdagpket(U, U, m, n, (1, p,)) 
    u5 = fp.docibrapdagpket(U, U, m, n, (p, n,)) 
    u6 = fp.docibrapdagpket(U, U, m, n, (n, p,)) 
    u7 = fp.docibrapdagpket(U, U, m, n, (p, m,)) 
    u8 = fp.docibrapdagpket(U, U, m, n, (m, p,)) 
    u9 = fp.docibrapdagpket(U, U, m, n, (q, 1,)) 
    u10 = fp.docibrapdagpket(U, U, m, n, (1, q,)) 
    u11 = fp.docibrapdagpket(U, U, m, n, (q, n,)) 
    u12 = fp.docibrapdagpket(U, U, m, n, (n, q,)) 
    u13 = fp.docibrapdagpket(U, U, m, n, (q, m,)) 
    u14 = fp.docibrapdagpket(U, U, m, n, (m, q,)) 
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

