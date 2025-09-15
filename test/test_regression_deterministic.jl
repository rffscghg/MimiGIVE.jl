module TestRegressionDeterministic

using MimiGIVE
include("utils.jl")

# label of validation data to compare AGAINST
validation_label = "current"

##------------------------------------------------------------------------------
## Validate Model Data
##------------------------------------------------------------------------------

savevars = [
    (compname=:DamageAggregator, varname=:total_damage),
    (compname=:DamageAggregator, varname=:total_damage_share),
    (compname=:DamageAggregator, varname=:total_damage_domestic),
    (compname=:DamageAggregator, varname=:cromar_mortality_damage),
    (compname=:DamageAggregator, varname=:agriculture_damage),
    (compname=:DamageAggregator, varname=:energy_damage),
    (compname=:DamageAggregator, varname=:cromar_mortality_damage_domestic),
    (compname=:DamageAggregator, varname=:agriculture_damage_domestic),
    (compname=:DamageAggregator, varname=:energy_damage_domestic),
    (compname=:global_netconsumption, varname=:net_consumption),
    (compname=:global_netconsumption, varname=:net_cpc),
    (compname=:global_netconsumption, varname=:global_gdp),
    (compname=:global_netconsumption, varname=:global_population),
    (compname=:temperature, varname=:T),
    (compname=:glaciers_small_icecaps, varname=:gsic_sea_level),
    (compname=:antarctic_icesheet, varname=:ais_sea_level),
    (compname=:greenland_icesheet, varname=:greenland_sea_level),
    (compname=:thermal_expansion, varname=:te_sea_level),
    (compname=:landwater_storage, varname=:lws_sea_level)
]

# default model
validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model")
m = MimiGIVE.get_model()
validate_model_data(m, savevars, validationdir)

# SSP245
validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model")
m = MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245")
validate_model_data(m, savevars, validationdir)

##------------------------------------------------------------------------------
## Validate SCC Data
##------------------------------------------------------------------------------
discount_rates = [
    # Constant discount rates
    (label="CR 1%", prtp=0.01, eta=0.0), (label="CR 2%", prtp=0.02, eta=0.0), (label="CR 2.5%", prtp=0.025, eta=0.0), (label="CR 3%", prtp=0.03, eta=0.0), (label="CR 5%", prtp=0.05, eta=0.0),
    # Some Ramsey discount rates
    (label="DICE2016", prtp=0.015, eta=1.45), (label="OtherRamsey", prtp=0.01, eta=1.)
]

for gas in [:CO2, :N2O, :CH4]
    # default model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
    validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model_SCC_2020")
    m = MimiGIVE.get_model()
    validate_scc_data(validationdir; m=m, year=2020, discount_rates=discount_rates, gas=gas)

    # SSP245 model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
    validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model_SCC_2020")
    m = MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245")
    validate_scc_data(validationdir; m=m, year=2020, discount_rates=discount_rates, gas=gas)
end


end # module
