using Pkg

# Instantiate environment
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Mimi
using MimiGIVE 
using DataFrames
using Query
using VegaLite

include("new_sector_damages.jl")
include("DamageAggregator_NewSectorDamages.jl")
include("main_model.jl")
include("mcs.jl")
include("scc.jl")

# Run the model
m = get_new_sector_model()
run(m)

# Explore the results in graphic form via the explorer and the plot functions
explore(m)
Mimi.plot(m, :NewSectorDamages, :damages)

# Examine results in tabular form, or plot them yourself with Vegalite
df = getdataframe(m, :NewSectorDamages, :damages) |> @filter(_.time >= 2020) |> DataFrame
df.time = string.(df.time)
df |> @vlplot(:line, x = "time:t", y = :damages, color = :country, width = 500, height = 500)

# Run a Monte Carlo Simulation
save_list = [
    (:DamageAggregator, :total_damage),
    (:DamageAggregator, :total_damage_share),
    (:DamageAggregator, :cromar_mortality_damage),
    (:DamageAggregator, :agriculture_damage),
    (:DamageAggregator, :energy_damage),
    (:DamageAggregator, :new_sector_damage),
    (:global_netconsumption, :net_consumption),
    (:global_netconsumption, :net_cpc),
    (:global_netconsumption, :global_gdp),
    (:global_netconsumption, :global_population),
    (:temperature, :T),
    (:glaciers_small_icecaps, :gsic_sea_level) ,
    (:antarctic_icesheet, :ais_sea_level),
    (:greenland_icesheet, :greenland_sea_level),
    (:thermal_expansion, :te_sea_level),
    (:landwater_storage, :lws_sea_level)
]

output_dir = joinpath(@__DIR__, "..", "output", "mcs", "MCS_main_output")
mkpath(output_dir)
results = run_new_sector_mcs(trials=10, save_list=save_list, output_dir=output_dir);

# Explore the results in graphic form via the explorer and the plot functions
explore(results)
Mimi.plot(results, :DamageAggregator, :new_sector_damage)

# Examine results in tabular form, or plot them yourself with Vegalite
getdataframe(results, :DamageAggregator, :new_sector_damage)

# Compute SCC
output_dir = joinpath(@__DIR__, "..", "output", "scc", "SCC_main_output")
mkpath(output_dir)
results = compute_new_sector_scc(year=2020, n=10, discount_rates = [(label = "DICE discount rate", prtp = 0.015, eta = 1.45)], output_dir=output_dir, save_list=save_list) # monte carlo
