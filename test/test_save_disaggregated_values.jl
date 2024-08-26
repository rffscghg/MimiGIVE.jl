module TestSaveDisaggregatedValues

using MimiGIVE
using Test
using Query
using DataFrames
using CSVFiles
using Mimi

atol = 1e-9
rtol = 1e-4

output_dir = mktempdir()
n = 3
m = MimiGIVE.get_model()
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45)]

results = MimiGIVE.compute_scc(m,
    year=2020,
    discount_rates=discount_rates,
    output_dir=output_dir,
    save_md=true,
    compute_sectoral_values=true,
    compute_disaggregated_values=true,
    compute_domestic_values=true,
    save_slr_damages=true,
    n=n,
    save_list=[(:Socioeconomic, :population),
        (:Socioeconomic, :gdp),
        (:DamageAggregator, :cromar_mortality_damage),
        (:DamageAggregator, :agriculture_damage),
        (:DamageAggregator, :energy_damage),
    ]
)

# Many of these tests check for internal consistency, meaning they check if the
# disaggregated results properly sum to the aggregated results, which are produced
# via the `save_list` or internal marginal damages calculations in the post trial 
# function. This will raise flags if either functionality changes without the other, 
# but checks via Figures are also encouraged.

##
## Sectoral Damages
##

# agriculture
df1 = DataFrame()
for region in dim_keys(m, :fund_regions)
    append!(df1,
        load(joinpath(output_dir, "results", "disaggregated_values", "damages_agriculture", "$(region).csv")) |>
        DataFrame |>
        i -> insertcols!(i, :region => Symbol(region))
    )
end
df1 = df1 |> @groupby({_.time, _.trialnum}) |> @map({key(_)..., damages = sum(_.damages)}) |> DataFrame
sort!(df1, [:trialnum, :time])
df2 = load(joinpath(output_dir, "results", "model_1", "DamageAggregator_agriculture_damage.csv")) |> @filter(_.time >= 2020) |> DataFrame
sort!(df2, [:trialnum, :time])
@test df1.damages ≈ df2.agriculture_damage

# mortality
df1 = DataFrame()
for country in dim_keys(m, :country)
    append!(df1,
        load(joinpath(output_dir, "results", "disaggregated_values", "damages_cromar_mortality", "$(country).csv")) |>
        DataFrame |>
        i -> insertcols!(i, :country => Symbol(country))
    )
end
df1 = df1 |> @groupby({_.time, _.trialnum}) |> @map({key(_)..., damages = sum(_.damages)}) |> DataFrame
sort!(df1, [:trialnum, :time])
df2 = load(joinpath(output_dir, "results", "model_1", "DamageAggregator_cromar_mortality_damage.csv")) |> @filter(_.time >= 2020) |> DataFrame
sort!(df2, [:trialnum, :time])
@test df1.damages ≈ df2.cromar_mortality_damage

# energy
df1 = DataFrame()
for country in dim_keys(m, :country)
    append!(df1,
        load(joinpath(output_dir, "results", "disaggregated_values", "damages_energy", "$(country).csv")) |>
        DataFrame |>
        i -> insertcols!(i, :country => Symbol(country))
    )
end
df1 = df1 |> @groupby({_.time, _.trialnum}) |> @map({key(_)..., damages = sum(_.damages)}) |> DataFrame
sort!(df1, [:trialnum, :time])
df2 = load(joinpath(output_dir, "results", "model_1", "DamageAggregator_energy_damage.csv")) |> @filter(_.time >= 2020) |> DataFrame
sort!(df2, [:trialnum, :time])
@test df1.damages ≈ df2.energy_damage

#slr
df1 = DataFrame()
for country in dim_keys(m, :country)
    filepath = joinpath(output_dir, "results", "disaggregated_values", "damages_slr", "$(country).csv")
    if isfile(filepath) # some countries not included
        append!(df1,
            load(filepath) |>
            DataFrame |>
            i -> insertcols!(i, :country => Symbol(country))
        )
    end
end
df1 = df1 |> @groupby({_.time, _.trialnum}) |> @map({key(_)..., damages = sum(_.damages)}) |> DataFrame
sort!(df1, [:trialnum, :time])
df2 = load(joinpath(output_dir, "results", "model_1", "slr_damages.csv")) |> @filter(_.time >= 2020) |> DataFrame
sort!(df2, [:trialnum, :time])
@test (df1.damages ./ 1e9) ≈ df2.slr_damages # convert disaggregated data into billions of USD

##
## Socioeconomics
##

# country_test = ["ABW", "CHI", "ZAF"] # random countries to test
gdp_savelist = load(joinpath(output_dir, "results", "model_1", "Socioeconomic_gdp.csv")) |>
               DataFrame |>
               @filter(_.time >= 2020) |>
               DataFrame

pop_savelist = load(joinpath(output_dir, "results", "model_1", "Socioeconomic_population.csv")) |>
               DataFrame |>
               @filter(_.time >= 2020) |>
               DataFrame

df_savelist = innerjoin(gdp_savelist, pop_savelist, on=[:time, :country, :trialnum])
insertcols!(df_savelist, :gdppc => df_savelist.gdp ./ df_savelist.population .* 1e3)
sort!(df_savelist, [:trialnum, :country, :time])

for country in Mimi.dim_keys(m, :country)
    disagg_values = load(joinpath(output_dir, "results", "disaggregated_values", "socioeconomics_country", "$country.csv")) |> DataFrame
    savelist_values = df_savelist |> @filter(_.country == country) |> DataFrame

    @test disagg_values.population ≈ savelist_values.population
    @test disagg_values.pc_gdp ≈ savelist_values.gdppc
end

##
## Marginal Damages 
## 

# compare domestic agriculture (FUND region) saved via original methods to the 
# disaggregated values added functionality

md_ag_domestic = DataFrame(results[:mds][(region=:domestic, sector=:agriculture)], Symbol.(2020:2300))
md_ag_domestic = md_ag_domestic |>
                 i -> insertcols!(i, 1, :trialnum => 1:3) |>
                      i -> stack(i, Not(:trialnum)) |>
                           i -> sort!(i, :trialnum) |>
                                DataFrame
md_usa = load(joinpath(output_dir, "results", "disaggregated_values", "mds_region_ag_only", "USA.csv")) |> DataFrame
@test md_usa.md ≈ md_ag_domestic.value

# compare domestic all other damages (USA + PRI) saved via original methods to the 
# disaggregated values added functionality

md_nonag_domestic = DataFrame(results[:mds][(region=:domestic, sector=:cromar_mortality)] .+ results[:mds][(region=:domestic, sector=:energy)] .+ results[:mds][(region=:domestic, sector=:slr)],
    Symbol.(2020:2300))
md_nonag_domestic = md_nonag_domestic |>
                    i -> insertcols!(i, 1, :trialnum => 1:3) |>
                         i -> stack(i, Not(:trialnum)) |>
                              i -> sort!(i, :trialnum) |>
                                   DataFrame

md_usa = load(joinpath(output_dir, "results", "disaggregated_values", "mds_country_no_ag", "USA.csv")) |> DataFrame
md_pri = load(joinpath(output_dir, "results", "disaggregated_values", "mds_country_no_ag", "PRI.csv")) |> DataFrame
md_usa_pri = copy(md_usa)
md_usa_pri.md = md_usa.md .+ md_pri.md
@test md_usa_pri.md ≈ md_nonag_domestic.value

end # module
