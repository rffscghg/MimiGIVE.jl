using MimiGIVE
using Random

include("utils.jl")

# This script saves a set of validation data in a post-fixed validation_label 
# subfolder of the validation_data folder.

# label of folder to be created
validation_label = "current"

##------------------------------------------------------------------------------
## Model Data
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
outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model")
isdir(outdir) || mkpath(outdir)
m = MimiGIVE.get_model()
save_model_data(m, savevars::Vector, outdir::String)

# SSP245
outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model")
isdir(outdir) || mkpath(outdir)
m = MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245")
save_model_data(m, savevars::Vector, outdir::String)

##------------------------------------------------------------------------------
## Compute SCC Data
##------------------------------------------------------------------------------
discount_rates = [
    # Constant discount rates
    (label="CR 1%", prtp=0.01, eta=0.0), (label="CR 2%", prtp=0.02, eta=0.0), (label="CR 2.5%", prtp=0.025, eta=0.0), (label="CR 3%", prtp=0.03, eta=0.0), (label="CR 5%", prtp=0.05, eta=0.0),
    # Some Ramsey discount rates
    (label="DICE2016", prtp=0.015, eta=1.45), (label="OtherRamsey", prtp=0.01, eta=1.)
]

# default model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model_SCC_2020")
isdir(outdir) || mkpath(outdir)

save_scc_data(outdir; m=MimiGIVE.get_model(), year=2020, discount_rates=discount_rates, gas=:CO2)
save_scc_data(outdir; m=MimiGIVE.get_model(), year=2020, discount_rates=discount_rates, gas=:CH4)
save_scc_data(outdir; m=MimiGIVE.get_model(), year=2020, discount_rates=discount_rates, gas=:N2O)

# SSP245, SC-CO2 and SC-CH4 and SC-N2O in year 2020
outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model_SCC_2020")
isdir(outdir) || mkpath(outdir)

save_scc_data(outdir; m=MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245"), year=2020, discount_rates=discount_rates, gas=:CO2)
save_scc_data(outdir; m=MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245"), year=2020, discount_rates=discount_rates, gas=:CH4)
save_scc_data(outdir; m=MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245"), year=2020, discount_rates=discount_rates, gas=:N2O)

##------------------------------------------------------------------------------
## Compute SCC MCS Data
##------------------------------------------------------------------------------

discount_rates = [
    # Constant discount rates
    (label="CR 1%", prtp=0.01, eta=0.0), (label="CR 2%", prtp=0.02, eta=0.0), (label="CR 2.5%", prtp=0.025, eta=0.0), (label="CR 3%", prtp=0.03, eta=0.0), (label="CR 5%", prtp=0.05, eta=0.0),
    # Some Ramsey discount rates
    (label="DICE2016", prtp=0.015, eta=1.45), (label="OtherRamsey", prtp=0.01, eta=1.)
]

save_list = [
    (:DamageAggregator, :total_damage),
    (:DamageAggregator, :total_damage_share),
    (:DamageAggregator, :total_damage_domestic),
    (:DamageAggregator, :cromar_mortality_damage),
    (:DamageAggregator, :agriculture_damage),
    (:DamageAggregator, :energy_damage),
    (:DamageAggregator, :cromar_mortality_damage_domestic),
    (:DamageAggregator, :agriculture_damage_domestic),
    (:DamageAggregator, :energy_damage_domestic),
    (:global_netconsumption, :net_consumption),
    (:global_netconsumption, :net_cpc),
    (:global_netconsumption, :global_gdp),
    (:global_netconsumption, :global_population),
    (:temperature, :T),
    (:glaciers_small_icecaps, :gsic_sea_level),
    (:antarctic_icesheet, :ais_sea_level),
    (:greenland_icesheet, :greenland_sea_level),
    (:thermal_expansion, :te_sea_level),
    (:landwater_storage, :lws_sea_level)
]

n = 3
seed = 999

# default model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
for gas in [:CO2, :CH4, :N2O]
    outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model_MCS_SCC_2020", "$gas")
    isdir(outdir) || mkpath(outdir)
    m = MimiGIVE.get_model()
    save_scc_mcs_data(seed, outdir, n; m=m, year=2020, discount_rates=discount_rates, gas=gas, save_list=save_list)
end

# SSP245, SC-CO2 and SC-CH4 and SC-N2O in year 2020
for gas in [:CO2, :CH4, :N2O]
    outdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model_MCS_SCC_2020", "$gas")
    isdir(outdir) || mkpath(outdir)
    m = MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245")
    save_scc_mcs_data(seed, outdir, n; m=m, year=2020, discount_rates=discount_rates, gas=gas, save_list=save_list)
end
