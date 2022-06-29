using Test, MimiGIVE

@testset "Mimi" begin

    @info("test_get_model.jl")
    @time include("test_get_model.jl")

    @info("test_compute_scc.jl")
    @time include("test_compute_scc.jl")
    
end
