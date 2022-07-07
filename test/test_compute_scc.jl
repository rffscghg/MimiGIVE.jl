# module TestComputeSCC

using MimiGIVE
using Test

import MimiGIVE: get_model, compute_scc

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
## API
##------------------------------------------------------------------------------

# deterministic - boilerplate API
scc = compute_scc(year = 2020)
scc_rff = compute_scc(get_model(; socioeconomics_source=:RFF, RFFSPsample=1000), year = 2020)
scc_ssp = compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year = 2020)

# deterministic - try non-default options, leave saving values for regression testing section
scc = compute_scc(year = 2025, last_year = 2200, prtp = 0.02, eta = 1.5, gas  = :CH4,
                CIAM_foresight = :limited, CIAM_GDPcap  = true, pulse_size = 10.)

# monte carlo - boilerplate API
drs = [(label = "label", prtp = 0.015, eta = 1.45)]
scc = compute_scc(year = 2020, n = 5, discount_rates = drs)
scc_rff = compute_scc(get_model(; socioeconomics_source=:RFF, RFFSPsample=1000), year = 2020, n = 5, discount_rates=drs)
scc_ssp = compute_scc(get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP126"), year = 2020, n = 5, discount_rates=drs)

# monte carlo - try non-default options, leave saving values for regression testing section
drs = [(label = "CR 1%", prtp = 0.01, eta = 0.0), 
        (label = "CR 2%", prtp = 0.02, eta = 0.0),
        (label = "CR 3%", prtp = 0.03, eta = 0.0)]
scc = compute_scc(year = 2025, last_year = 2200, discount_rates = drs, gas  = :CH4,
        CIAM_foresight = :limited, CIAM_GDPcap  = true, pulse_size = 10., n = 5)

##------------------------------------------------------------------------------
## keyword arguments and values
##------------------------------------------------------------------------------

# year and last_year
@test compute_scc(year = 2020) < compute_scc(year = 2025) < compute_scc(year = 2030)
@test compute_scc(year = 2020) > compute_scc(year = 2020; last_year=2200)

# discount rates parameters eta and prtp
@test compute_scc(year = 2020, prtp = 0.01, eta = 0.0) > compute_scc(year = 2020, prtp = 0.02, eta = 0.0) > compute_scc(year = 2020, prtp = 0.03, eta = 0.0) 
drs = [(label = "CR 1%", prtp = 0.01, eta = 0.0), 
        (label = "CR 2%", prtp = 0.02, eta = 0.0),
        (label = "CR 3%", prtp = 0.03, eta = 0.0)]
sccs = compute_scc(year = 2020; discount_rates = drs)
@test sccs[(dr_label = "CR 1%", prtp = 0.01, eta = 0.0)] > sccs[(dr_label = "CR 2%", prtp = 0.02, eta = 0.0)] > sccs[(dr_label = "CR 3%", prtp = 0.03, eta = 0.0)]

# gas
@test compute_scc(year = 2020, gas = :CO2) < compute_scc(year = 2020, gas = :CH4) < compute_scc(year = 2020, gas = :N2O)

# pulse size
scc_1 = compute_scc(year = 2020, pulse_size = 1.)
scc_1_1 = compute_scc(year = 2020, pulse_size = 1.1)

# end # module