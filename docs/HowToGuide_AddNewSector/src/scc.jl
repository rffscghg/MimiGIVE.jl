using Mimi
using MimiGIVE
using Distributions
using Dates
using Query

"""
Compute the SC of a gas in USD \$2005
"""
function compute_modified_scc(
    m::Model=get_modified_model();
    year::Union{Int,Nothing}=nothing,
    last_year::Int=MimiGIVE._model_years[end],
    prtp::Union{Float64,Nothing}=0.015,
    eta::Union{Float64,Nothing}=1.45,
    discount_rates=nothing,
    certainty_equivalent=false,
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Vector{Int},Nothing}=nothing,
    rffsp_sampling::Symbol=:random,
    rffsp_sampling_ids::Union{Vector{Int},Nothing}=nothing,
    n=0,
    gas::Symbol=:CO2,
    save_list::Vector=[],
    output_dir::Union{String,Nothing}=nothing,
    save_md::Bool=false,
    save_cpc::Bool=false,
    save_slr_damages::Bool=false,
    compute_sectoral_values::Bool=false,
    compute_domestic_values::Bool=false,
    CIAM_foresight::Symbol=:perfect,
    CIAM_GDPcap::Bool=false,
    post_mcs_creation_function=nothing,
    pulse_size::Float64=1.
)

    hfc_list = [:HFC23, :HFC32, :HFC43_10, :HFC125, :HFC134a, :HFC143a, :HFC227ea, :HFC245fa]
    gases_list = [:CO2, :CH4, :N2O, hfc_list...]

    m = deepcopy(m) # in the case that an `m` was provided, be careful that we don't modify the original

    year === nothing ? error("Must specify an emission year. Try `compute_scc(m, year=2020)`.") : nothing
    !(last_year in MimiGIVE._model_years) ? error("Invalid value of $last_year for last_year. last_year must be within the model's time index $(MimiGIVE._model_years).") : nothing
    !(year in MimiGIVE._model_years) ? error("Cannot compute the scc for year $year, year must be within the model's time index $(MimiGIVE._model_years).") : nothing
    !(gas in gases_list) ? error("Invalid value of $gas for gas, gas must be one of $(gases_list).") : nothing
    n > 0 && certainty_equivalent && !save_cpc && error("certainty_equivalent=true also requires save_cpc=true")

    mm = MimiGIVE.get_marginal_model(m; year=year, gas=gas, pulse_size=pulse_size)

    if n == 0
        return MimiGIVE._compute_scc(mm,
            year=year,
            last_year=last_year,
            prtp=prtp,
            eta=eta,
            discount_rates=discount_rates,
            gas=gas,
            domestic=compute_domestic_values,
            CIAM_foresight=CIAM_foresight,
            CIAM_GDPcap=CIAM_GDPcap,
            pulse_size=pulse_size
        )
    else
        isnothing(discount_rates) ? error("To run the Monte Carlo compute_scc function (n != 0), please use the `discount_rates` argument.") : nothing

        # Set up output directories
        output_dir = output_dir === nothing ? joinpath(@__DIR__, "../output/mcs-SC/", "MCS $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$n") : output_dir
        isdir("$output_dir/results") || mkpath("$output_dir/results")

        return _compute_modified_scc_mcs(mm,
            n,
            year=year,
            last_year=last_year,
            discount_rates=discount_rates,
            certainty_equivalent=certainty_equivalent,
            fair_parameter_set=fair_parameter_set,
            fair_parameter_set_ids=fair_parameter_set_ids,
            rffsp_sampling=rffsp_sampling,
            rffsp_sampling_ids=rffsp_sampling_ids,
            gas=gas,
            save_list=save_list,
            output_dir=output_dir,
            save_md=save_md,
            save_cpc=save_cpc,
            save_slr_damages=save_slr_damages,
            compute_sectoral_values=compute_sectoral_values,
            compute_domestic_values=compute_domestic_values,
            CIAM_foresight=CIAM_foresight,
            CIAM_GDPcap=CIAM_GDPcap,
            post_mcs_creation_function=post_mcs_creation_function,
            pulse_size=pulse_size
        )
    end
end

# Post trial function to to after each trial within the MCS
function modified_post_trial_func(mcs::SimulationInstance, trialnum::Int, ntimesteps::Int, tup)

    # Unpack the payload object 
    scc_values, intermediate_ce_scc_values, md_values, cpc_values, slr_damages, year, last_year, discount_rates, gas, ciam_base, ciam_modified, segment_fingerprints, options = Mimi.payload2(mcs)

    # Compute some useful indices
    year_index = findfirst(isequal(year), MimiGIVE._model_years)
    last_year_index = findfirst(isequal(last_year), MimiGIVE._model_years)

    # Access the models
    base, marginal = mcs.models
    damages_base = base[:DamageAggregator, :total_damage]
    damages_marginal = marginal[:DamageAggregator, :total_damage]

    if options.compute_domestic_values
        damages_base_domestic = base[:DamageAggregator, :total_damage_domestic]
        damages_marginal_domestic = marginal[:DamageAggregator, :total_damage_domestic]
    end

    # Compute marginal damages
    # Units Note:
    #   main_mds and non-ciam sectoral damages: we explicitly need to handle both pulse size and molecular mass so we use gas_units_multiplier
    #   slr_mds: within the _compute_ciam_marginal_damages function we handle both pulse size and molecular mass
    gas_units_multiplier = MimiGIVE.scc_gas_molecular_conversions[gas] ./ (MimiGIVE.scc_gas_pulse_size_conversions[gas] .* options.pulse_size)
    include_slr = base[:DamageAggregator, :include_slr]

    if include_slr
        ciam_mds = MimiGIVE._compute_ciam_marginal_damages(base, marginal, gas, ciam_base, ciam_modified, segment_fingerprints; CIAM_foresight=options.CIAM_foresight, CIAM_GDPcap=options.CIAM_GDPcap, pulse_size=options.pulse_size) # NamedTuple with globe and domestic
        # zero out the CIAM marginal damages from start year (2020) through emissions
        # year - they will be non-zero due to foresight but saved marginal damages
        # should be zeroed out pre-emissions year
        ciam_mds.globe[1:year_index] .= 0.
        ciam_mds.domestic[1:year_index] .= 0.
    end

    main_mds = (damages_marginal .- damages_base) .* gas_units_multiplier
    slr_mds = include_slr ? ciam_mds.globe : fill(0., length(MimiGIVE._model_years))
    total_mds = main_mds .+ slr_mds

    if options.compute_domestic_values
        main_mds_domestic = (damages_marginal_domestic .- damages_base_domestic) .* gas_units_multiplier
        slr_mds_domestic = include_slr ? ciam_mds.domestic : fill(0., length(MimiGIVE._model_years))
        total_mds_domestic = main_mds_domestic .+ slr_mds_domestic
    end

    if options.compute_sectoral_values

        cromar_mortality_mds = (marginal[:DamageAggregator, :cromar_mortality_damage] .- base[:DamageAggregator, :cromar_mortality_damage]) .* gas_units_multiplier
        agriculture_mds = (marginal[:DamageAggregator, :agriculture_damage] .- base[:DamageAggregator, :agriculture_damage]) .* gas_units_multiplier
        energy_mds = (marginal[:DamageAggregator, :energy_damage] .- base[:DamageAggregator, :energy_damage]) .* gas_units_multiplier
        new_sector_mds = (marginal[:DamageAggregator, :new_sector_damage] .- base[:DamageAggregator, :new_sector_damage]) .* gas_units_multiplier

        if options.compute_domestic_values
            cromar_mortality_mds_domestic = (marginal[:DamageAggregator, :cromar_mortality_damage_domestic] .- base[:DamageAggregator, :cromar_mortality_damage_domestic]) .* gas_units_multiplier
            agriculture_mds_domestic = (marginal[:DamageAggregator, :agriculture_damage_domestic] .- base[:DamageAggregator, :agriculture_damage_domestic]) .* gas_units_multiplier
            energy_mds_domestic = (marginal[:DamageAggregator, :energy_damage_domestic] .- base[:DamageAggregator, :energy_damage_domestic]) .* gas_units_multiplier
            new_sector_mds_domestic = (marginal[:DamageAggregator, :new_sector_damage_domestic] .- base[:DamageAggregator, :new_sector_damage_domestic]) .* gas_units_multiplier
        end
    end

    # Save marginal damages
    if options.save_md

        # global
        md_values[(region=:globe, sector=:total)][trialnum, :] = total_mds[MimiGIVE._damages_idxs]
        if options.compute_sectoral_values
            md_values[(region=:globe, sector=:cromar_mortality)][trialnum, :] = cromar_mortality_mds[MimiGIVE._damages_idxs]
            md_values[(region=:globe, sector=:agriculture)][trialnum, :] = agriculture_mds[MimiGIVE._damages_idxs]
            md_values[(region=:globe, sector=:energy)][trialnum, :] = energy_mds[MimiGIVE._damages_idxs]
            md_values[(region=:globe, sector=:slr)][trialnum, :] = slr_mds[MimiGIVE._damages_idxs]
            md_values[(region=:globe, sector=:new_sector)][trialnum, :] = new_sector_mds[MimiGIVE._damages_idxs]
        end

        # domestic
        if options.compute_domestic_values
            md_values[(region=:domestic, sector=:total)][trialnum, :] = total_mds_domestic[MimiGIVE._damages_idxs]
            if options.compute_sectoral_values
                md_values[(region=:domestic, sector=:cromar_mortality)][trialnum, :] = cromar_mortality_mds_domestic[MimiGIVE._damages_idxs]
                md_values[(region=:domestic, sector=:agriculture)][trialnum, :] = agriculture_mds_domestic[MimiGIVE._damages_idxs]
                md_values[(region=:domestic, sector=:energy)][trialnum, :] = energy_mds_domestic[MimiGIVE._damages_idxs]
                md_values[(region=:domestic, sector=:slr)][trialnum, :] = slr_mds_domestic[MimiGIVE._damages_idxs]
                md_values[(region=:domestic, sector=:new_sector)][trialnum, :] = new_sector_mds_domestic[MimiGIVE._damages_idxs]
            end
        end
    end

    # Save slr damages
    if options.save_slr_damages

        # get a dummy ciam model to be sure to accurately assign segment names to 
        # segment level damages
        m = get_modified_model()
        m_ciam, ~ = MimiGIVE.get_ciam(m)

        if include_slr

            # global
            slr_damages[:base][trialnum, :] = ciam_mds.damages_base[MimiGIVE._damages_idxs]
            slr_damages[:modified][trialnum, :] = ciam_mds.damages_modified[MimiGIVE._damages_idxs]
            slr_damages[:base_lim_cnt][trialnum, :, :] = ciam_mds.base_lim_cnt
            slr_damages[:modified_lim_cnt][trialnum, :, :] = ciam_mds.modified_lim_cnt
            slr_damages[:base_segments_2100][trialnum, :] = ciam_mds.damages_base_segments_2100

            # domestic - these Dictionary entries will only exist if we are computing
            # domestic values
            if options.compute_domestic_values
                slr_damages[:base_domestic][trialnum, :] = ciam_mds.damages_base_domestic[MimiGIVE._damages_idxs]
                slr_damages[:modified_domestic][trialnum, :] = ciam_mds.damages_modified_domestic[MimiGIVE._damages_idxs]
            end

        else

            # global
            slr_damages[:base][trialnum, :] .= 0.
            slr_damages[:modified][trialnum, :] .= 0.
            slr_damages[:base_lim_cnt][trialnum, :, :] .= 0.
            slr_damages[:modified_lim_cnt][trialnum, :, :] .= 0.
            slr_damages[:base_segments_2100][trialnum, :] .= 0.

            # domestic - these Dictionary entries will only exist if we are computing
            # domestic values
            if options.compute_domestic_values
                slr_damages[:base_domestic][trialnum, :] .= 0.
                slr_damages[:modified_domestic][trialnum, :] .= 0.
            end
        end
    end

    # Get per capita consumption
    # We don't care about units here because we are only going to use ratios
    cpc = base[:global_netconsumption, :net_cpc]

    # Save per capita consumption
    if options.save_cpc
        cpc_values[(region=:globe, sector=:total)][trialnum, :] = cpc[MimiGIVE._damages_idxs]
    end

    # Calculate the SCC for each discount rate
    for dr in discount_rates
        df = [((cpc[year_index] / cpc[i])^dr.eta * 1 / (1 + dr.prtp)^(t - year) for (i, t) in enumerate(MimiGIVE._model_years) if year <= t <= last_year)...]
        if options.certainty_equivalent
            df_ce = [((1. / cpc[i])^dr.eta * 1 / (1 + dr.prtp)^(t - year) for (i, t) in enumerate(MimiGIVE._model_years) if year <= t <= last_year)...] # only used if optionas.certainty_equivalent=true
        end

        # totals (sector=:total)
        scc = sum(df .* total_mds[year_index:last_year_index])
        scc_values[(region=:globe, sector=:total, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc
        if options.certainty_equivalent
            intermediate_ce_scc = sum(df_ce .* total_mds[year_index:last_year_index])
            intermediate_ce_scc_values[(region=:globe, sector=:total, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
        end

        # domestic totals (sector=:total)
        if options.compute_domestic_values
            scc = sum(df .* total_mds_domestic[year_index:last_year_index])
            scc_values[(region=:domestic, sector=:total, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            if options.certainty_equivalent
                intermediate_ce_scc = sum(df_ce .* total_mds_domestic[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:domestic, sector=:total, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
            end
        end

        # sectoral
        if options.compute_sectoral_values
            scc = sum(df .* cromar_mortality_mds[year_index:last_year_index])
            scc_values[(region=:globe, sector=:cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            scc = sum(df .* agriculture_mds[year_index:last_year_index])
            scc_values[(region=:globe, sector=:agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            scc = sum(df .* energy_mds[year_index:last_year_index])
            scc_values[(region=:globe, sector=:energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            scc = sum(df .* slr_mds[year_index:last_year_index])
            scc_values[(region=:globe, sector=:slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            scc = sum(df .* new_sector_mds[year_index:last_year_index])
            scc_values[(region=:globe, sector=:new_sector, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

            if options.certainty_equivalent
                intermediate_ce_scc = sum(df_ce .* cromar_mortality_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                intermediate_ce_scc = sum(df_ce .* agriculture_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                intermediate_ce_scc = sum(df_ce .* energy_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                intermediate_ce_scc = sum(df_ce .* slr_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                intermediate_ce_scc = sum(df_ce .* new_sector_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:new_sector, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
            end

            # sectoral domestic (region=:domestic)
            if options.compute_domestic_values

                scc = sum(df .* cromar_mortality_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector=:cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

                scc = sum(df .* agriculture_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector=:agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

                scc = sum(df .* energy_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector=:energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

                scc = sum(df .* slr_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector=:slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

                scc = sum(df .* new_sector_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector=:new_sector, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc

                if options.certainty_equivalent
                    intermediate_ce_scc = sum(df_ce .* cromar_mortality_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector=:cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                    intermediate_ce_scc = sum(df_ce .* agriculture_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector=:agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                    intermediate_ce_scc = sum(df_ce .* energy_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector=:energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                    intermediate_ce_scc = sum(df_ce .* slr_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector=:slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                    intermediate_ce_scc = sum(df_ce .* new_sector_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector=:new_sector, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc

                end
            end
        end
    end
end

# Internal function to compute the SCC in a Monte Carlo Simulation
function _compute_modified_scc_mcs(mm::MarginalModel,
    n;
    year::Int,
    last_year::Int,
    discount_rates,
    certainty_equivalent::Bool,
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Vector{Int},Nothing}=nothing,
    rffsp_sampling::Symbol=:random,
    rffsp_sampling_ids::Union{Vector{Int},Nothing}=nothing,
    gas::Symbol,
    save_list::Vector,
    output_dir::String,
    save_md::Bool,
    save_cpc::Bool,
    save_slr_damages::Bool,
    compute_sectoral_values::Bool,
    compute_domestic_values::Bool,
    CIAM_foresight::Symbol,
    CIAM_GDPcap::Bool,
    post_mcs_creation_function,
    pulse_size::Float64
)

    models = [mm.base, mm.modified]

    socioeconomics_module = MimiGIVE._get_module_name(mm.base, :Socioeconomic)
    if socioeconomics_module == :MimiSSPs
        socioeconomics_source = :SSP
    elseif socioeconomics_module == :MimiRFFSPs
        socioeconomics_source = :RFF
    end

    mcs = get_modified_mcs(n;
        socioeconomics_source=socioeconomics_source,
        mcs_years=MimiGIVE._model_years,
        fair_parameter_set=fair_parameter_set,
        fair_parameter_set_ids=fair_parameter_set_ids,
        rffsp_sampling=rffsp_sampling,
        rffsp_sampling_ids=rffsp_sampling_ids,
        save_list=save_list
    )

    if post_mcs_creation_function !== nothing
        post_mcs_creation_function(mcs)
    end

    regions = compute_domestic_values ? [:globe, :domestic] : [:globe]
    sectors = compute_sectoral_values ? [:total, :cromar_mortality, :agriculture, :energy, :slr, :new_sector] : [:total]

    scc_values = Dict((region=r, sector=s, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta) => Vector{Float64}(undef, n) for dr in discount_rates, r in regions, s in sectors)
    intermediate_ce_scc_values = certainty_equivalent ? Dict((region=r, sector=s, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta) => Vector{Float64}(undef, n) for dr in discount_rates, r in regions, s in sectors) : nothing
    md_values = save_md ? Dict((region=r, sector=s) => Array{Float64}(undef, n, length(MimiGIVE._damages_years)) for r in regions, s in sectors) : nothing
    cpc_values = save_cpc ? Dict((region=r, sector=s) => Array{Float64}(undef, n, length(MimiGIVE._damages_years)) for r in [:globe], s in [:total]) : nothing # just global and total for now

    if save_slr_damages

        # global
        slr_damages = Dict(
            :base => Array{Float64}(undef, n, length(MimiGIVE._damages_years)),
            :modified => Array{Float64}(undef, n, length(MimiGIVE._damages_years)),
            :base_lim_cnt => Array{Float64}(undef, n, length(MimiGIVE._damages_years), 145), # 145 CIAM countries
            :modified_lim_cnt => Array{Float64}(undef, n, length(MimiGIVE._damages_years), 145), # 145 CIAM countries
            :base_segments_2100 => Array{Float64}(undef, n, 11835) # 11,835 segments
        )

        # domestic
        # optionally add arrays to hold the domestic base and modified damages
        if compute_domestic_values
            slr_damages[:base_domestic] = Array{Float64}(undef, n, length(MimiGIVE._damages_years))
            slr_damages[:modified_domestic] = Array{Float64}(undef, n, length(MimiGIVE._damages_years))
        end
    else
        slr_damages = nothing
    end

    ciam_base, segment_fingerprints = MimiGIVE.get_ciam(mm.base)
    ciam_modified, _ = MimiGIVE.get_ciam(mm.base)

    ciam_base = Mimi.build(ciam_base)
    ciam_modified = Mimi.build(ciam_modified)

    # set some computation options
    options = (
        compute_sectoral_values=compute_sectoral_values,
        compute_domestic_values=compute_domestic_values,
        save_md=save_md,
        save_cpc=save_cpc,
        save_slr_damages=save_slr_damages,
        CIAM_foresight=CIAM_foresight,
        CIAM_GDPcap=CIAM_GDPcap,
        certainty_equivalent=certainty_equivalent,
        pulse_size=pulse_size
    )

    payload = [scc_values, intermediate_ce_scc_values, md_values, cpc_values, slr_damages, year, last_year, discount_rates, gas, ciam_base, ciam_modified, segment_fingerprints, options]

    Mimi.set_payload2!(mcs, payload)

    # Run all model years even if taking a shorter last_year - running unnecessary 
    # timesteps but simplifies accumulation     
    sim_results = run(mcs,
        models,
        n,
        post_trial_func=modified_post_trial_func,
        results_in_memory=false,
        results_output_dir="$output_dir/results"
    )

    # unpack the payload object
    scc_values, intermediate_ce_scc_values, md_values, cpc_values, slr_damages, year, last_year, discount_rates, gas, ciam_base, ciam_modified, segment_fingerprints, options = Mimi.payload2(sim_results)

    # Write out the slr damages to disk in the same place that variables from the save_list would be written out
    if save_slr_damages
        isdir("$output_dir/results/model_1") || mkpath("$output_dir/results/model_1")
        isdir("$output_dir/results/model_2") || mkpath("$output_dir/results/model_2")

        # global 
        df = DataFrame(slr_damages[:base], :auto) |>
             i -> rename!(i, Symbol.(MimiGIVE._damages_years)) |>
                  i -> insertcols!(i, 1, :trial => 1:n) |>
                       i -> stack(i, Not(:trial)) |>
                            i -> rename!(i, [:trial, :time, :slr_damages]) |>
                                 save("$output_dir/results/model_1/slr_damages.csv")

        df = DataFrame(slr_damages[:modified], :auto) |>
             i -> rename!(i, Symbol.(MimiGIVE._damages_years)) |>
                  i -> insertcols!(i, 1, :trial => 1:n) |>
                       i -> stack(i, Not(:trial)) |>
                            i -> rename!(i, [:trial, :time, :slr_damages]) |>
                                 save("$output_dir/results/model_2/slr_damages.csv")

        segments = Symbol.(dim_keys(ciam_base, :segments))
        df = DataFrame(slr_damages[:base_segments_2100], :auto) |>
             i -> rename!(i, segments) |>
                  i -> insertcols!(i, 1, :trial => 1:n) |>
                       i -> stack(i, Not(:trial)) |>
                            i -> rename!(i, [:trial, :segment, :slr_damages_2100]) |>
                                 save("$output_dir/results/model_1/slr_damages_2100_by_segment.csv")

        # domestic 
        if compute_domestic_values
            df = DataFrame(slr_damages[:base_domestic], :auto) |>
                 i -> rename!(i, Symbol.(MimiGIVE._damages_years)) |>
                      i -> insertcols!(i, 1, :trial => 1:n) |>
                           i -> stack(i, Not(:trial)) |>
                                i -> rename!(i, [:trial, :time, :slr_damages_domestic]) |>
                                     save("$output_dir/results/model_1/slr_damages_domestic.csv")

            df = DataFrame(slr_damages[:modified_domestic], :auto) |>
                 i -> rename!(i, Symbol.(MimiGIVE._damages_years)) |>
                      i -> insertcols!(i, 1, :trial => 1:n) |>
                           i -> stack(i, Not(:trial)) |>
                                i -> rename!(i, [:trial, :time, :slr_damages_domestic]) |>
                                     save("$output_dir/results/model_2/slr_damages_domestic.csv")
        end

        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))

        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))
        df = DataFrame(:trial => [], :time => [], :country => [], :capped_flag => [])
        for trial in 1:n # loop over trials
            trial_df = DataFrame(slr_damages[:base_lim_cnt][trial, :, :], :auto) |>
                       i -> rename!(i, ciam_country_names) |>
                            i -> insertcols!(i, 1, :time => MimiGIVE._damages_years) |>
                                 i -> stack(i, Not(:time)) |>
                                      i -> insertcols!(i, 1, :trial => fill(trial, length(MimiGIVE._damages_years) * 145)) |>
                                           i -> rename!(i, [:trial, :time, :country, :capped_flag]) |>
                                                i -> @filter(i, _.capped_flag == 1) |>
                                                     DataFrame
            append!(df, trial_df)
        end
        df |> save("$output_dir/results/slr_damages_base_lim_counts.csv")

        df = DataFrame(:trial => [], :time => [], :country => [], :capped_flag => [])
        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))
        for trial in 1:n # loop over trials
            trial_df = DataFrame(slr_damages[:modified_lim_cnt][trial, :, :], :auto) |>
                       i -> rename!(i, ciam_country_names) |>
                            i -> insertcols!(i, 1, :time => MimiGIVE._damages_years) |>
                                 i -> stack(i, Not(:time)) |>
                                      i -> insertcols!(i, 1, :trial => fill(trial, length(MimiGIVE._damages_years) * 145)) |>
                                           i -> rename!(i, [:trial, :time, :country, :capped_flag]) |>
                                                i -> @filter(i, _.capped_flag == 1) |>
                                                     DataFrame
            append!(df, trial_df)
        end
        df |> save("$output_dir/results/slr_damages_modified_lim_counts.csv")
    end

    expected_mu_in_year_of_emission = Dict()

    if certainty_equivalent
        year_index = findfirst(isequal(year), MimiGIVE._damages_years)
        # In this case the normalization from utils to $ hasn't happened in the post trial function
        # and instead we now do this here, based on expected per capita consumption in the year
        # of the marginal emission pulse
        cpc_in_year_of_emission = view(cpc_values[(region=:globe, sector=:total)], :, year_index)

        for k in keys(scc_values)
            expected_mu_in_year_of_emission[k] = mean(1 ./ (cpc_in_year_of_emission .^ k.eta))
        end
    end

    # Construct the returned result object
    result = Dict()

    # add an :scc dictionary, where key value pairs (k,v) are NamedTuples with keys(prtp, eta, region, sector) => values are 281 element vectors (2020:2300)
    result[:scc] = Dict()
    for (k, v) in scc_values
        if certainty_equivalent
            result[:scc][k] = (
                expected_scc=mean(v),
                se_expected_scc=std(v) / sqrt(n),
                ce_scc=mean(intermediate_ce_scc_values[k]) ./ expected_mu_in_year_of_emission[k],
                ce_sccs=intermediate_ce_scc_values[k] ./ expected_mu_in_year_of_emission[k],
                sccs=v,
            )
        else
            result[:scc][k] = (
                expected_scc=mean(v),
                se_expected_scc=std(v) / sqrt(n),
                sccs=v
            )
        end
    end

    # add a :mds dictionary, where key value pairs (k,v) are NamedTuples with keys(region, sector) => values are (n x 281 (2020:2300)) matrices
    if save_md
        result[:mds] = Dict()
        for (k, v) in md_values
            result[:mds][k] = v
        end
    end

    # add a :cpc dictionary, where key value pairs (k,v) are NamedTuples with keys(region, sector) => values are (n x 281 (2020:2300)) matrices
    if save_cpc
        result[:cpc] = Dict()
        for (k, v) in cpc_values
            result[:cpc][k] = v
        end
    end

    return result
end
