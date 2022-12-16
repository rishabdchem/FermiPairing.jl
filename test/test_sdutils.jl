using FermiPairing 
const fp = FermiPairing 
using Test

@testset "Indexing" begin

    # Initialize
    m = rand(4:40)
    n = rand(1:m-1)
    i = rand(1:n)
    a = rand(n+1:m)

    # Tests
    HF = zeros(Int, m)
    HF[1:n] = ones(Int, n)
    @test fp.sdgen(HF, n, (a,), (i,)) == fp.sdindex(fp.phgen(m, n, (a,), (i,)), n)

end

