module TestSaveDisaggregatedValues

using MimiGIVE
using Test
using Query
using DataFrames
using CSVFiles

atol = 1e-9

output_dir = mktempdir()
# output_dir = "/Users/lisarennels/.julia/dev/MimiGIVE/test/tmp"
n = 3
m = MimiGIVE.get_model()
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45)]

results = MimiGIVE.compute_scc(m, 
                                year=2020, 
                                discount_rates = discount_rates,
                                output_dir = output_dir,
                                save_md = true,
                                compute_sectoral_values = true,
                                compute_disaggregated_values = true,
                                compute_domestic_values = true,
                                n = n
                            )

# test damages 

# test socioeconomics

##
## Test Marginal Damages 
## 

# compare domestic agriculture (FUND region) saved via original methods to the 
# disaggregated values added functionality

md_ag_domestic = DataFrame(results[:mds][(region = :domestic, sector = :agriculture)], Symbol.(2020:2300))
md_ag_domestic = md_ag_domestic |>
                    i -> insertcols!(i, 1, :trial => 1:3) |>
                    i -> stack(i, Not(:trial)) |>
                    i -> sort!(i, :trial) |>
                    DataFrame
md_usa = load(joinpath(output_dir, "results", "disaggregated_values", "mds_region_ag_only", "USA.csv")) |> DataFrame
@test md_usa.md ≈ md_ag_domestic.value atol = atol

# compare domestic all other damages (USA + PRI) saved via original methods to the 
# disaggregated values added functionality

md_nonag_domestic = DataFrame(results[:mds][(region = :domestic, sector = :cromar_mortality)] .+ results[:mds][(region = :domestic, sector = :energy)] .+ results[:mds][(region = :domestic, sector = :slr)], 
                                Symbol.(2020:2300))
md_nonag_domestic = md_nonag_domestic |>
                i -> insertcols!(i, 1, :trial => 1:3) |>
                i -> stack(i, Not(:trial)) |>
                i -> sort!(i, :trial) |>
                DataFrame

md_usa = load(joinpath(output_dir, "results", "disaggregated_values", "mds_country_no_ag", "USA.csv")) |> DataFrame
md_pri = load(joinpath(output_dir, "results", "disaggregated_values", "mds_country_no_ag", "PRI.csv")) |> DataFrame
md_usa_pri = copy(md_usa)
md_usa_pri.md = md_usa.md .+ md_pri.md
@test md_usa_pri.md ≈ md_nonag_domestic.value atol = atol

end # module