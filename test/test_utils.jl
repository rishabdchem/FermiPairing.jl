using FermiPairing 
const fp = FermiPairing 
using LinearAlgebra 
using Test

@testset "ESP" begin

    # Initialize
    m = rand(1:100)
    V = randn(m)

    # ESP
    e1 = getesp(V, 1) 
    e2 = getesp(V, m) 

    # Tests
    t1 = 0.0
    t2 = 1.0
    for i = 1:m
        t1 += V[i]    
        t2 *= V[i]    
    end
    @test e1 == t1
    @test e2 == t2

end


@testset "NOCI" begin

    # Initialize
    m = rand(1:100)
    M = Symmetric(randn(m, m))

    # Tests
    e1, V1, zm = fp.nocigs(M) 
    D, V = eigen(M)
    e2, V2, = D[1], V[:, 1] 
    @test isapprox(e1, e2, atol = 1e-12) 
    @test isapprox(V1, V2, atol = 1e-12) 

end


