using Test, MimiGIVE
ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

@testset "Mimi" begin

    @info("test_get_model.jl")
    @time include("test_get_model.jl")

    @info("test_compute_scc.jl")
    @time include("test_compute_scc.jl")

    @info("test_regression_deterministic.jl")
    @time include("test_regression_deterministic.jl")

    if VERSION == v"1.10" # random number generator not always stable between versions
        @info("test_regression_mcs.jl")
        @time include("test_regression_mcs.jl")
    end 

    @info("test_disaggregated_values.jl")
    @time include("test_save_disaggregated_values.jl")
    
end

