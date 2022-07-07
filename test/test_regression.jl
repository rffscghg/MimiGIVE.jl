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
validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_date", "default_model")
m = MimiGIVE.get_model()
validate_model_data(m, savevars, validationdir)

# SSP245
validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_date", "SSP245_model")
m = MimiGIVE.get_model(; socioeconomics_source = :SSP, SSP_scenario = "SSP245")
validate_model_data(m, savevars, validationdir)

##------------------------------------------------------------------------------
## Compute SCC Data
##------------------------------------------------------------------------------
discount_rates = [
                    # Constant discount rates
                    (label = "CR 1%", prtp = 0.01, eta = 0.0), (label = "CR 2%", prtp = 0.02, eta = 0.0), (label = "CR 2.5%", prtp = 0.025, eta = 0.0), (label = "CR 3%", prtp = 0.03, eta = 0.0), (label = "CR 5%", prtp = 0.05, eta = 0.0),
                    # Ramsey discount rates calibrated a la NPP
                    (label = "1.5%", prtp = exp(9.149606e-05) - 1, eta = 1.016010e+00), (label = "2.0%", prtp = exp(0.001972641) - 1, eta = 1.244458999), (label = "2.5%", prtp = exp(0.004618784) - 1, eta = 1.421158088), (label = "3.0%", prtp = exp(0.007702711) - 1, eta = 1.567899391),
                    # Other Ramsey discount rates
                    (label = "DICE2016", prtp = 0.015, eta = 1.45)
                ]

# default model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_date", "default_model_SCC_2020")
m = MimiGIVE.get_model()
validate_scc_data(validationdir; m = m, year = 2020, discount_rates = discount_rates, gas = :CO2)

# SSP245 model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
validationdir = joinpath(@__DIR__, "validation_data","validation_data_$validation_date", "SSP245_model_SCC_2020")
m = MimiGIVE.get_model(; socioeconomics_source = :SSP, SSP_scenario = "SSP245")
validate_scc_data(validationdir; m = m, year = 2020, discount_rates = discount_rates, gas = :CO2)

##------------------------------------------------------------------------------
## Compute SCC MCS Data
##------------------------------------------------------------------------------

end # module
