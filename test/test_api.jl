module TestAPI

using MimiGIVE
using Test

# This testing module will test that the public facing API (1) runs without error
# and (2) produces outputs of the expected dimensions and types

##------------------------------------------------------------------------------
# get_model
##------------------------------------------------------------------------------

# function get_model(; Agriculture_gtap::String = "midDF",
#     socioeconomics_source::Symbol = :RFF,
#     SSP_scenario::Union{Nothing, String} = nothing,       
#     RFFSPsample::Union{Nothing, Int} = nothing,
#     Agriculture_floor_on_damages::Bool = true,
#     Agriculture_ceiling_on_benefits::Bool = false,
#     vsl::Symbol= :epa
# )

m = MimiGIVE.get_model()
run(m)

##------------------------------------------------------------------------------
# compute_scc
##------------------------------------------------------------------------------

function compute_scc(m::Model=get_model(); 
    year::Union{Int, Nothing} = nothing, 
    last_year::Int = _model_years[end], 
    prtp::Union{Float64,Nothing} = 0.015, 
    eta::Union{Float64,Nothing}=1.45,
    discount_rates=nothing,
    certainty_equivalent=false,
    fair_parameter_set::Symbol = :random,
    fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
    rffsp_sampling::Symbol = :random,
    rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
    n=0,
    gas::Symbol = :CO2,
    save_list::Vector = [],
    output_dir::Union{String, Nothing} = nothing,
    save_md::Bool = false,
    save_cpc::Bool = false,
    save_slr_damages::Bool = false,
    compute_sectoral_values::Bool = false,
    compute_domestic_values::Bool = false,
    CIAM_foresight::Symbol = :perfect,
    CIAM_GDPcap::Bool = false,
    post_mcs_creation_function=nothing,
    pulse_size::Float64=1.
)

end # module
