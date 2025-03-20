using FileIO, CSVFiles, DataFrames, Query

# Helper functions for streaming disaggregated data within the compute_scc funciton

function _stream_disagg_damages(m::Mimi.Model, output_dir::String, trialnum::Int, streams::Dict)
    # println("Streaming out trialnum $trialnum ...")	
    cromar_mortality_damages = getdataframe(m, :DamageAggregator, :damage_cromar_mortality) |>
                               @filter(_.time > 2019) |>
                               @rename(:damage_cromar_mortality => :damages) |>
                               DataFrame |>
                               i -> insertcols!(i, 1, :trialnum => trialnum)

    energy_damages = getdataframe(m, :DamageAggregator, :damage_energy) |>
                     @filter(_.time > 2019) |>
                     @mutate(damage_energy = _.damage_energy * 1e9) |> # billions of USD to USD	
                     @rename(:damage_energy => :damages, :energy_countries => :country) |>
                     DataFrame |>
                     i -> insertcols!(i, 1, :trialnum => trialnum)

    ag_damages = getdataframe(m, :DamageAggregator, :damage_ag) |>
                 @filter(_.time > 2019) |>
                 @mutate(damage_ag = _.damage_ag * 1e9) |> # billions of USD to USD	
                 @rename(:damage_ag => :damages, :fund_regions => :region) |>
                 DataFrame |>
                 i -> insertcols!(i, 1, :trialnum => trialnum)

    for country in unique(cromar_mortality_damages.country)
        filename = joinpath("$output_dir/results/disaggregated_values/damages_cromar_mortality/$(country).csv")
        trial_df = cromar_mortality_damages |> @filter(_.country == country) |> @select(:trialnum, :time, :damages) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end

    for country in unique(energy_damages.country)
        filename = joinpath("$output_dir/results/disaggregated_values/damages_energy/$(country).csv")
        trial_df = energy_damages |> @filter(_.country == country) |> @select(:trialnum, :time, :damages) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end

    for region in unique(ag_damages.region)
        filename = joinpath("$output_dir/results/disaggregated_values/damages_agriculture/$(region).csv")
        trial_df = ag_damages |> @filter(_.region == region) |> @select(:trialnum, :time, :damages) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end
end

function _stream_disagg_socioeconomics(m::Mimi.Model, output_dir::String, trialnum::Int, streams::Dict)
    # println("Streaming out trialnum $trialnum ...")	

    country_pop = getdataframe(m, :Socioeconomic, :population) |> @filter(_.time > 2019) |> DataFrame # millions	
    country_pc_gdp = getdataframe(m, :PerCapitaGDP, :pc_gdp) |> @filter(_.time > 2019) |> DataFrame  # USD per capita	
    country_data = innerjoin(country_pop, country_pc_gdp, on=[:time, :country]) |> i -> insertcols!(i, 1, :trialnum => trialnum)

    for country in unique(country_data.country)
        filename = joinpath("$output_dir/results/disaggregated_values/socioeconomics_country/$(country).csv")
        trial_df = country_data |> @filter(_.country == country) |> @select(:trialnum, :time, :population, :pc_gdp) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end

    region_pop = getdataframe(m, :Agriculture, :population) |> @filter(_.time > 2019) |> DataFrame # millions	
    region_gdp = getdataframe(m, :Agriculture, :income) |> @filter(_.time > 2019) |> DataFrame # billions USD	
    region_data = innerjoin(region_pop, region_gdp, on=[:time, :fund_regions])
    region_data = insertcols!(region_data, :pc_gdp => (region_data.income ./ region_data.population) * 1e3) |>
                  @select(:time, :fund_regions, :population, :pc_gdp) |>
                  DataFrame |>
                  i -> insertcols!(i, 1, :trialnum => trialnum)

    for region in unique(region_data.fund_regions)
        filename = joinpath("$output_dir/results/disaggregated_values/socioeconomics_region/$(region).csv")
        trial_df = region_data |> @filter(_.fund_regions == region) |> @select(:trialnum, :time, :population, :pc_gdp) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end
end

# note we pass a ModelInstance because ciam_base and ciam_modified are instances	
function _stream_disagg_damages_slr(m::Mimi.ModelInstance, data::Array, output_dir::String, trialnum::Int, streams::Dict)
    slr_damages = DataFrame(data, dim_keys(m, :ciam_country)) |>
                  i -> insertcols!(i, 1, :time => _damages_years) |>
                       i -> stack(i, Not(:time)) |>
                            @filter(_.time > 2019) |>
                            @rename(:variable => :country, :value => :damages) |>
                            @mutate(damages = _.damages * 1e9) |> # billions of USD to USD	
                            DataFrame |>
                            i -> insertcols!(i, 1, :trialnum => trialnum)


    for country in unique(slr_damages.country)
        filename = joinpath("$output_dir/results/disaggregated_values/damages_slr/$(country).csv")
        trial_df = slr_damages |> @filter(_.country == country) |> @select(:trialnum, :time, :damages) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end
end

function _stream_disagg_md(m_base::Mimi.Model, m_modified::Mimi.Model, ciam_base::Union{Nothing,Mimi.ModelInstance}, md_ciam::Union{Nothing,Array},
    output_dir::String, trialnum::Int, streams::Dict; gas_units_multiplier::Float64)

    # get marginal damages in USD $2005 and be sure to adjust for # adjust for the (1) molecular mass and (2) pulse size, as well as billions of USD to USD for ag and energy	
    md_cromar_mortality = (view(m_modified[:DamageAggregator, :damage_cromar_mortality], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_cromar_mortality], _damages_idxs, :)) .* gas_units_multiplier
    md_energy = (view(m_modified[:DamageAggregator, :damage_energy], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_energy], _damages_idxs, :)) .* 1e9 .* gas_units_multiplier
    md_ag = (view(m_modified[:DamageAggregator, :damage_ag], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_ag], _damages_idxs, :)) .* 1e9 .* gas_units_multiplier

    # save agriculture	
    md_ag_df = DataFrame(md_ag, dim_keys(m_base, :fund_regions)) |>
               i -> insertcols!(i, 1, :time => _damages_years) |>
                    i -> stack(i, Not(:time)) |>
                         @rename(:variable => :region, :value => :md) |>
                         DataFrame |>
                         i -> insertcols!(i, 1, :trialnum => trialnum)

    for region in unique(md_ag_df.region)
        filename = joinpath("$output_dir/results/disaggregated_values/mds_region_ag_only/$(region).csv")
        trial_df = md_ag_df |> @filter(_.region == region) |> @select(:trialnum, :time, :md) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end

    # save country level mds	

    # aggregate ciam marginal damages	
    md_ciam_all_countries = fill(0., size(md_cromar_mortality))
    if !isnothing(ciam_base)
        country_idxs = indexin(dim_keys(ciam_base, :ciam_country), dim_keys(m_base, :country))
        md_ciam_all_countries[:, country_idxs] = md_ciam[_damages_idxs, :]
    end

    md_country_df = DataFrame(md_cromar_mortality .+ md_energy .+ md_ciam_all_countries, dim_keys(m_base, :country)) |>
                    i -> insertcols!(i, 1, :time => _damages_years) |>
                         i -> stack(i, Not(:time)) |>
                              @rename(:variable => :country, :value => :md) |>
                              DataFrame |>
                              i -> insertcols!(i, 1, :trialnum => trialnum)

    for country in unique(md_country_df.country)
        filename = joinpath("$output_dir/results/disaggregated_values/mds_country_no_ag/$(country).csv")
        trial_df = md_country_df |> @filter(_.country == country) |> @select(:trialnum, :time, :md) |> DataFrame
        if haskey(streams, filename)
            write(streams[filename], trial_df)
        else
            streams[filename] = savestreaming(filename, trial_df)
        end
    end
end
