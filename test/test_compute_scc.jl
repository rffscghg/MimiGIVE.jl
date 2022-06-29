module TestComputeSCC

using MimiGIVE
using Test

import MimiGIVE: get_model, compute_scc

# This module will test that the public facing API (1) runs without error
# (2) picks up relevant keyword arguments and (2) produces outputs of the 
# expected dimensions and types

##------------------------------------------------------------------------------
## compute_scc
##------------------------------------------------------------------------------

# function compute_scc(m::Model=get_model(); 
#     year::Union{Int, Nothing} = nothing, 
#     last_year::Int = _model_years[end], 
#     prtp::Union{Float64,Nothing} = 0.015, 
#     eta::Union{Float64,Nothing}=1.45,
#     discount_rates=nothing,
#     certainty_equivalent=false,
#     fair_parameter_set::Symbol = :random,
#     fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
#     rffsp_sampling::Symbol = :random,
#     rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
#     n=0,
#     gas::Symbol = :CO2,
#     save_list::Vector = [],
#     output_dir::Union{String, Nothing} = nothing,
#     save_md::Bool = false,
#     save_cpc::Bool = false,
#     save_slr_damages::Bool = false,
#     compute_sectoral_values::Bool = false,
#     compute_domestic_values::Bool = false,
#     CIAM_foresight::Symbol = :perfect,
#     CIAM_GDPcap::Bool = false,
#     post_mcs_creation_function=nothing,
#     pulse_size::Float64=1.
# )

##------------------------------------------------------------------------------
## basic API
##------------------------------------------------------------------------------

# deterministic
scc = compute_scc(year = 2020)
scc_rff = compute_scc(get_model(; socioeconomics_source=:RFF, RFFSPsample=1000), year = 2020)
scc_ssp = compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year = 2020)

# monte carlo
drs = [(label = "label", prtp = 0.015, eta = 1.45)]
scc = compute_scc(year = 2020, n = 5, discount_rates = drs)
scc_rff = compute_scc(get_model(; socioeconomics_source=:RFF, RFFSPsample=1000), year = 2020, n = 5, discount_rates=drs)
scc_ssp = compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year = 2020, n = 5, discount_rates=drs)

##------------------------------------------------------------------------------
## keyword arguments and values
##------------------------------------------------------------------------------

end # module