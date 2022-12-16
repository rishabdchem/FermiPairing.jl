using Test

@testset "FermiPairing tests" begin
    include("test_utils.jl")
    include("test_sdutils.jl")
    include("test_flexible.jl")
    include("test_sdci.jl")
    include("test_agputils.jl")
    include("test_agp.jl")
    include("test_lcagp.jl")
end
