using MimiGIVE
include("utils.jl")

# This script saves a set of validation data for the give date in the validation_data
# folder.  These data are used for regression testing by test_regression.jl.

curr_date = "07_07_2022"

##------------------------------------------------------------------------------
## Model Data
##------------------------------------------------------------------------------

savevars = [
    (compname = :DamageAggregator, varname = :total_damage),
    (compname = :DamageAggregator, varname = :total_damage_share),
    (compname = :DamageAggregator, varname = :total_damage_domestic),
    (compname = :DamageAggregator, varname = :cromar_mortality_damage),
    (compname = :DamageAggregator, varname = :agriculture_damage),
    (compname = :DamageAggregator, varname = :energy_damage),
    (compname = :DamageAggregator, varname = :cromar_mortality_damage_domestic),
    (compname = :DamageAggregator, varname = :agriculture_damage_domestic),
    (compname = :DamageAggregator, varname = :energy_damage_domestic),
    (compname = :global_netconsumption, varname = :net_consumption),
    (compname = :global_netconsumption, varname = :net_cpc),
    (compname = :global_netconsumption, varname = :global_gdp),
    (compname = :global_netconsumption, varname = :global_population),   
]

# default model
outdir = joinpath(@__DIR__, "validation_data", "default_model", "validation_data_$curr_date")
isdir(outdir) || mkpath(outdir)
m = MimiGIVE.get_model()
save_model_data(m, savevars::Vector, outdir::String)

# SSP245
outdir = joinpath(@__DIR__, "validation_data", "SSP245_model", "validation_data_$curr_date")
isdir(outdir) || mkpath(outdir)
m = MimiGIVE.get_model(; socioeconomics_source = :SSP, SSP_scenario = "SSP245")
save_model_data(m, savevars::Vector, outdir::String)

##------------------------------------------------------------------------------
## Compute SCC Data
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## Compute SCC MCS Data
##------------------------------------------------------------------------------
