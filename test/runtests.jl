using Test, MimiGIVE

@testset "Mimi" begin

    @info("test_api.jl")
    @time include("test_api.jl")
    
end
