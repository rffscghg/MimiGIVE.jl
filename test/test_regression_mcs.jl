module TestRegressionMCS

using MimiGIVE
include("utils.jl")

# label of validation data to compare AGAINST
validation_label = "current"

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
for gas in [:CO2, :N2O, :CH4]
    validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "default_model_MCS_SCC_2020", "$gas")
    m = MimiGIVE.get_model()
    validate_scc_mcs_data(seed, validationdir, n;
        m=m,
        year=2020,
        discount_rates=discount_rates,
        gas=gas,
        save_list=save_list,
        save_md=true,
        save_cpc=true,
        save_slr_damages=true,
        compute_sectoral_values=true,
        compute_domestic_values=true,
    )
end

# SSP245 model, SC-CO2 and SC-CH4 and SC-N2O in year 2020
for gas in [:CO2, :N2O, :CH4]
    validationdir = joinpath(@__DIR__, "validation_data", "validation_data_$validation_label", "SSP245_model_MCS_SCC_2020", "$gas")
    m = MimiGIVE.get_model(; socioeconomics_source=:SSP, SSP_scenario="SSP245")
    validate_scc_mcs_data(seed, validationdir, n;
        m=m,
        year=2020,
        discount_rates=discount_rates,
        gas=gas,
        save_list=save_list,
        save_md=true,
        save_cpc=true,
        save_slr_damages=true,
        compute_sectoral_values=true,
        compute_domestic_values=true,
    )
end

end # module
