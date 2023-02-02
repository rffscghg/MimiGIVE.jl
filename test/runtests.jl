using Test, MimiGIVE
ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

@testset "Mimi" begin

    @info("test_get_model.jl")
    @time include("test_get_model.jl")

    @info("test_compute_scc.jl")
    @time include("test_compute_scc.jl")

    @info("test_regression.jl")
    @time include("test_regression.jl")
    
end
