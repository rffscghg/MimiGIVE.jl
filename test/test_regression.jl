module TestRegression

using MimiGIVE
include("utils.jl")

# date of validation data
validation_date = "07_07_2022"

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
validationdir = joinpath(@__DIR__, "validation_data", "default_model", "validation_data_$validation_date")
m = MimiGIVE.get_model()
validate_model_data(m, savevars, validationdir)

# SSP245
validationdir = joinpath(@__DIR__, "validation_data", "SSP245_model", "validation_data_$validation_date")
m = MimiGIVE.get_model(; socioeconomics_source = :SSP, SSP_scenario = "SSP245")
validate_model_data(m, savevars, validationdir)

##------------------------------------------------------------------------------
## Compute SCC Data
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## Compute SCC MCS Data
##------------------------------------------------------------------------------

end # module
