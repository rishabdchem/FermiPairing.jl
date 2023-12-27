using Test

@testset "FermiPairing tests" begin
    include("test_utils.jl")
    include("test_sdutils.jl")
    include("test_flexible.jl")
    include("test_sdci.jl")
    include("test_agputils.jl")
    include("test_agp.jl")
    include("test_lcagp.jl")
    include("test_optlcagp.jl")
    include("test_apig.jl")
    include("test_btsutils.jl")
    include("test_btsoverlap.jl")
    include("test_singlebts.jl")
    include("test_btstranrdm.jl")
    include("test_lcbts.jl")
    include("test_optlcbts.jl")
end
