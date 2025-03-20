using Pkg

# Instantiate environment
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Mimi
using MimiGIVE
using DataFrames
using Query
using VegaLite

include("main_model.jl")
include("main_mcs.jl")
include("scc.jl")

# Run the model
m = get_modified_model()
run(m)

# Explore the results in graphic form via the explorer and the plot functions
explore(m)
Mimi.plot(m, :NewSectorDamages, :damages)

# Examine results in tabular form, or plot them yourself with Vegalite
df = getdataframe(m, :NewSectorDamages, :damages) |> @filter(_.time >= 2020) |> DataFrame
df.time = string.(df.time)
df |> @vlplot(:line, x = "time:t", y = :damages, color = :country, width = 500, height = 500)

# Run a Monte Carlo Simulation - NOTE be careful, at a high value of n any
# country-disaggregated variable will save to an extremely large file
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
    (:glaciers_small_icecaps, :gsic_sea_level),
    (:antarctic_icesheet, :ais_sea_level),
    (:greenland_icesheet, :greenland_sea_level),
    (:thermal_expansion, :te_sea_level),
    (:landwater_storage, :lws_sea_level)
]

output_dir = joinpath(@__DIR__, "..", "output", "mcs", "MCS_main_output")
mkpath(output_dir)
results = run_modified_mcs(trials=10, save_list=save_list, output_dir=output_dir);

# Explore the results in graphic form via the explorer and the plot functions
explore(results)
Mimi.plot(results, :DamageAggregator, :new_sector_damage)

# Examine results in tabular form, or plot them yourself with Vegalite
getdataframe(results, :DamageAggregator, :new_sector_damage)

# Compute SCC
output_dir = joinpath(@__DIR__, "..", "output", "scc", "SCC_main_output")
mkpath(output_dir)
results = compute_modified_scc(
    year=2020,
    n=10,
    discount_rates=[
        (label="DICE discount rate", prtp=0.015, eta=1.45),
        (label="2.0%", prtp=exp(0.001972641) - 1, eta=1.244458999)
    ],
    output_dir=output_dir,
    save_list=save_list,
    compute_sectoral_values=true,
    save_md=true
);

# Access the computed SCC values
scc_df = DataFrame(:region => [], :sector => [], :discount_rate_label => [], :expected_scc => [], :se_expected_scc => [], :ew => [], :ew_norm_region => [])
results_scc = results[:scc] # results is a dictionary, :scc is a key to this dictionary

for (k, v) in results_scc
    # results_scc is a dictionary, we iterate over it's keys (k) and values (v)
    # --- the keys are each a NamedTuple with elements region, sector, dr_label, prtp, and eta
    # --- the values are each a Named Tuple with elements expected_scc, se_expected_scc, and sccs (a vector of the sccs)
    append!(scc_df, DataFrame(
        :region => k.region,
        :sector => k.sector,
        :ew => k.ew,
        :ew_norm_region => k.ew_norm_region,
        :discount_rate_label => k.dr_label,
        :expected_scc => v[:expected_scc],
        :se_expected_scc => v[:se_expected_scc]
    )
    )
end
show(scc_df)

# Access the marginal damages (undiscounted)

# marginal damages for the global region for sector new_sector
mds = results[:mds][((region=:globe, sector=:new_sector))]
mds_df = DataFrame(mds, :auto)
rename!(mds_df, Symbol.(2020:2300))
insertcols!(mds_df, 1, :trial => 1:10)
mds_df = stack(mds_df, Not(:trial))
rename!(mds_df, [:trial, :time, :value])

mds_df |> @vlplot(:line, x = "time:t", y = :value, color = "trial:n")
