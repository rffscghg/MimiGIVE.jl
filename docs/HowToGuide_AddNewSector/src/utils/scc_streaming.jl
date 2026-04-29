using FileIO, CSVFiles, DataFrames, Query

# Helper functions for streaming disaggregated data within the compute_scc funciton

function _modified_stream_disagg_damages(m::Mimi.Model, output_dir::String, trialnum::Int, streams::Dict)
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

    new_sector_damages = getdataframe(m, :DamageAggregator, :damage_new_sector) |>
                         @filter(_.time > 2019) |>
                         @mutate(damage_new_sector = _.damage_new_sector * 1e9) |> # billions of USD to USD	
                         @rename(:damage_new_sector => :damages, :new_sector_countries => :country) |>
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

    for country in unique(new_sector_damages.country)
        filename = joinpath("$output_dir/results/disaggregated_values/damages_new_sector/$(country).csv")
        trial_df = new_sector_damages |> @filter(_.country == country) |> @select(:trialnum, :time, :damages) |> DataFrame
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

function _modified_stream_disagg_md(m_base::Mimi.Model, m_modified::Mimi.Model, ciam_base::Union{Nothing,Mimi.ModelInstance}, md_ciam::Union{Nothing,Array},
    output_dir::String, trialnum::Int, streams::Dict; gas_units_multiplier::Float64)

    # get marginal damages in USD $2005 and be sure to adjust for # adjust for the (1) molecular mass and (2) pulse size, as well as billions of USD to USD for ag and energy	
    md_cromar_mortality = (view(m_modified[:DamageAggregator, :damage_cromar_mortality], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_cromar_mortality], _damages_idxs, :)) .* gas_units_multiplier
    md_energy = (view(m_modified[:DamageAggregator, :damage_energy], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_energy], _damages_idxs, :)) .* 1e9 .* gas_units_multiplier
    md_new_sector = (view(m_modified[:DamageAggregator, :damage_new_sector], _damages_idxs, :) .- view(m_base[:DamageAggregator, :damage_new_sector], _damages_idxs, :)) .* 1e9 .* gas_units_multiplier
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

    md_country_df = DataFrame(md_cromar_mortality .+ md_energy .+ md_new_sector .+ md_ciam_all_countries, dim_keys(m_base, :country)) |>
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
