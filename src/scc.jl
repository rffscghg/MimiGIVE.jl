using Dates, CSVFiles, DataFrames

const _model_years = collect(1750:2300)
const _damages_years = collect(2020:2300)
const _damages_idxs = indexin(_damages_years, _model_years)

const scc_gas_molecular_conversions = Dict(:CO2 => 12/44, # C to CO2
                                            :N2O => 28/44, # N2 to N2O,
                                            :CH4 => 1., # CH4 to CH4
                                            :HFC23 => 1., # HFC23 to HFC23
                                            :HFC32 => 1., # HFC32 to HFC32
                                            :HFC43_10 => 1., # HFC43_10 to HFC43_10
                                            :HFC125 => 1., # HFC125 to HFC125
                                            :HFC134a => 1., # HFC134a to HFC134a
                                            :HFC143a => 1., # HFC143a to HFC143a
                                            :HFC227ea => 1., # HFC227ea to HFC227ea
                                            :HFC245fa => 1.) # HFC245fa to HFC245fa

const scc_gas_pulse_size_conversions = Dict(:CO2 => 1e9, # Gt to t
                                        :N2O => 1e6, # Mt to t
                                        :CH4 => 1e6, # Mt to t
                                        :HFC23 => 1e3, # kt to t
                                        :HFC32 => 1e3, # kt to t
                                        :HFC43_10 => 1e3, # kt to t
                                        :HFC125 => 1e3, # kt to t
                                        :HFC134a => 1e3, # kt to t
                                        :HFC143a => 1e3, # kt to t
                                        :HFC227ea => 1e3, # kt to t
                                        :HFC245fa => 1e3) # kt to t
"""
    compute_scc(m::Model=get_model(); 
            year::Union{Int, Nothing} = nothing,
            last_year::Int = _model_years[end],
            prtp::Union{Float64,Nothing} = 0.015,
            eta::Union{Float64,Nothing}=1.45,
            discount_rates=nothing,
            certainty_equivalent=false,
            fair_parameter_set::Symbol = :random,
            fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
            rffsp_sampling::Symbol = :random,
            rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
            n=0,
            gas::Symbol = :CO2,
            save_list::Vector = [],
            output_dir::Union{String, Nothing} = nothing,
            save_md::Bool = false,
            save_cpc::Bool = false,
            save_slr_damages::Bool = false,
            compute_sectoral_values::Bool = false,
            compute_domestic_values::Bool = false,
            CIAM_foresight::Symbol = :perfect,
            CIAM_GDPcap::Bool = false,
            post_mcs_creation_function=nothing,
            pulse_size::Float64=1.
        )

Compute the SC of a gas for the GIVE in USD \$2005

- `m` (default get_model()) - If no model is provided, the default model from MimiGIVE.get_model() is used. 
- `prtp` (default 0.015) and `eta` (1.45) - Ramsey discounting parameterization
- `discount_rates` (default nothing) - a vector of Named Tuples ie. [(prpt = 0.03., eta = 1.45), (prtp = 0.015, eta = 1.45)] - required if running n > 1
- `certainty_equivalent` (default false) - whether to compute the certainty equivalent or expected SCC
- `fair_parameter_set` (default :random) - :random means FAIR mcs samples will be 
chosen randomly from the provided sets, while :deterministic means they will be 
based on the provided vector of to `fair_parameter_set_ids` keyword argument. 
- `fair_parameter_set_ids` - (default nothing) - if `fair_parameter_set` is set 
to :deterministic, this `n` element vector provides the fair parameter set ids 
that will be run, otherwise it is set to `nothing` and ignored.
- `rffsp_sampling` (default :random) - which sampling strategy to use for the RFF 
SPs, :random means RFF SPs will be chosen randomly, while :deterministic means they 
will be based on the provided vector of to `rffsp_sampling_ids` keyword argument. 
- `rffsp_sampling_ids` - (default nothing) - if `rffsp_sampling` is set to :deterministic, 
this `n` element vector provides the RFF SP ids that will be run, otherwise it is 
set to `nothing` and ignored.
- `n` (default 0) - If `n` is 0, the deterministic version will be run, otherwise, a monte carlo simulation will be run. 
- `gas` (default :CO2) - the gas for which to compute the SC, options are :CO2, :CH4, and :N2O. 
- `save_list` (default []) - which parameters and varaibles to save for each trial,
entered as a vector of Tuples (:component_name, :variable_name)
- `output_dir` (default constructed folder name) - folder to hold results 
- `save_md` (default is false) - save and return the marginal damages from a monte carlo simulation
- `save_cpc` (default is false) - save and return the per capita consumption from a monte carlo simulation
- `save_slr_damages`(default is false) - save global sea level rise damages from CIAM to disk
- `compute_sectoral_values` (default is false) - compute and return sectoral values as well as total
- `compute_domestic_values` (default is false) - compute and return domestic values in addition to global
- CIAM_foresight(default is :perfect) - Use limited foresight (:limited) or perfect foresight (:perfect) for MimiCIAM cost calculations
- CIAM_GDPcap (default is false) - Limit SLR damages to country-level annual GDP
- `pulse_size` (default 1.) - This determines the size of the additional pulse of emissions. Default of `1.` implies the standard pulse size 
of 1Gt of C for CO2, 1Mt of CH4, and 1Mt of N2O. 

"""
function compute_scc(m::Model=get_model(); 
            year::Union{Int, Nothing} = nothing, 
            last_year::Int = _model_years[end], 
            prtp::Union{Float64,Nothing} = 0.015, 
            eta::Union{Float64,Nothing}=1.45,
            discount_rates=nothing,
            certainty_equivalent=false,
            fair_parameter_set::Symbol = :random,
            fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
            rffsp_sampling::Symbol = :random,
            rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
            n=0,
            gas::Symbol = :CO2,
            save_list::Vector = [],
            output_dir::Union{String, Nothing} = nothing,
            save_md::Bool = false,
            save_cpc::Bool = false,
            save_slr_damages::Bool = false,
            compute_sectoral_values::Bool = false,
            compute_domestic_values::Bool = false,
            CIAM_foresight::Symbol = :perfect,
            CIAM_GDPcap::Bool = false,
            post_mcs_creation_function=nothing,
            pulse_size::Float64=1.
        )

    hfc_list = [:HFC23, :HFC32, :HFC43_10, :HFC125, :HFC134a, :HFC143a, :HFC227ea, :HFC245fa]
    gases_list = [:CO2, :CH4, :N2O, hfc_list ...]

    m = deepcopy(m) # in the case that an `m` was provided, be careful that we don't modify the original

    year === nothing ? error("Must specify an emission year. Try `compute_scc(m, year=2020)`.") : nothing
    !(last_year in _model_years) ? error("Invalid value of $last_year for last_year. last_year must be within the model's time index $_model_years.") : nothing
    !(year in _model_years) ? error("Cannot compute the scc for year $year, year must be within the model's time index $_model_years.") : nothing
    !(gas in gases_list) ? error("Invalid value of $gas for gas, gas must be one of $(gases_list).") : nothing
    n>0 && certainty_equivalent && !save_cpc && error("certainty_equivalent=true also requires save_cpc=true")
    
    mm = get_marginal_model(m; year = year, gas = gas, pulse_size = pulse_size)

    if n==0
        return _compute_scc(mm, 
                            year=year, 
                            last_year=last_year, 
                            prtp=prtp, 
                            eta=eta, 
                            discount_rates=discount_rates, 
                            gas=gas, domestic=compute_domestic_values, 
                            CIAM_foresight=CIAM_foresight,
                            CIAM_GDPcap=CIAM_GDPcap,
                            pulse_size=pulse_size
                        )
    else
        isnothing(discount_rates) ? error("To run the Monte Carlo compute_scc function (n != 0), please use the `discount_rates` argument.") : nothing
        
        # Set up output directories
        output_dir = output_dir === nothing ? joinpath(@__DIR__, "../output/mcs-SC/", "MCS $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$n") : output_dir
        isdir("$output_dir/results") || mkpath("$output_dir/results")

        return _compute_scc_mcs(mm, 
                                n, 
                                year=year, 
                                last_year=last_year, 
                                discount_rates=discount_rates, 
                                certainty_equivalent=certainty_equivalent,
                                fair_parameter_set = fair_parameter_set,
                                fair_parameter_set_ids = fair_parameter_set_ids,
                                rffsp_sampling = rffsp_sampling,
                                rffsp_sampling_ids = rffsp_sampling_ids,
                                gas = gas, 
                                save_list = save_list, 
                                output_dir = output_dir,
                                save_md = save_md,
                                save_cpc = save_cpc,
                                save_slr_damages = save_slr_damages,
                                compute_sectoral_values = compute_sectoral_values,
                                compute_domestic_values = compute_domestic_values,
                                CIAM_foresight = CIAM_foresight,
                                CIAM_GDPcap = CIAM_GDPcap,
                                post_mcs_creation_function = post_mcs_creation_function,
                                pulse_size=pulse_size
                            )
    end
end

# helper function for computing SCC from a MarginalModel, not to be exported or advertised to users
function _compute_scc(mm::MarginalModel; 
                        year::Int, 
                        last_year::Int, 
                        prtp, 
                        eta, 
                        discount_rates, 
                        gas::Symbol,
                        domestic::Bool,
                        CIAM_foresight::Symbol,
                        CIAM_GDPcap::Bool,
                        pulse_size::Float64
                    )
                    
    # Run all model years even if taking a shorter last_year - running unnecessary 
    # timesteps but simplifies accumulation             
    run(mm)

    # at this point create identical copies ciam_base and ciam_modified, they will 
    # be updated in _compute_ciam_marginal_damages with update_ciam!
    ciam_base, segment_fingerprints = get_ciam(mm.base)
    ciam_modified, _ = get_ciam(mm.base) 

    ciam_base = Mimi.build(ciam_base)
    ciam_modified = Mimi.build(ciam_modified)

    # Units Note:
    #   main_marginal_damages: the marginal model will handle pulse size, we handle molecular mass conversion explicilty
    #   ciam_marginal_damages: within the _compute_ciam_marginal_damages function we handle both pulse size and molecular mass
    if domestic
        main_marginal_damages = mm[:DamageAggregator, :total_damage_domestic] .* scc_gas_molecular_conversions[gas] 
        ciam_marginal_damages = mm.base[:DamageAggregator, :include_slr] ? _compute_ciam_marginal_damages(mm.base, mm.modified, gas, ciam_base, ciam_modified, segment_fingerprints; CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,  pulse_size=pulse_size).domestic : fill(0., length(_model_years)) 
    else
        main_marginal_damages = mm[:DamageAggregator, :total_damage] .* scc_gas_molecular_conversions[gas] 
        ciam_marginal_damages = mm.base[:DamageAggregator, :include_slr] ? _compute_ciam_marginal_damages(mm.base, mm.modified, gas, ciam_base, ciam_modified, segment_fingerprints; CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap, pulse_size=pulse_size).globe : fill(0., length(_model_years)) 
    end

    marginal_damages = main_marginal_damages .+ ciam_marginal_damages
    
    # We don't care about units here because we are only going to use ratios
    cpc = mm.base[:global_netconsumption, :net_cpc]

    year_index = findfirst(isequal(year), _model_years)
    last_year_index = findfirst(isequal(last_year), _model_years)

    if discount_rates!==nothing
        sccs = Dict{NamedTuple{(:dr_label, :prtp,:eta),Tuple{Any, Float64,Float64}}, Float64}()

        for dr in discount_rates
            df = [((cpc[year_index]/cpc[i])^dr.eta * 1/(1+dr.prtp)^(t-year) for (i,t) in enumerate(_model_years) if year<=t<=last_year)...]
            scc = sum(df .* marginal_damages[year_index:last_year_index])

            sccs[(dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)] = scc
        end

        return sccs
    else
        df = [((cpc[year_index]/cpc[i])^eta * 1/(1+prtp)^(t-year) for (i,t) in enumerate(_model_years) if year<=t<=last_year)...]
        scc = sum(df .* marginal_damages[year_index:last_year_index])

        return scc
    end
end

function post_trial_func(mcs::SimulationInstance, trialnum::Int, ntimesteps::Int, tup)

    # Unpack the payload object 
    scc_values, intermediate_ce_scc_values, md_values, cpc_values, slr_damages, year, last_year, discount_rates, gas, ciam_base, ciam_modified, segment_fingerprints, options = Mimi.payload2(mcs)

    # Compute some useful indices
    year_index = findfirst(isequal(year), _model_years)
    last_year_index = findfirst(isequal(last_year), _model_years)

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
    gas_units_multiplier = scc_gas_molecular_conversions[gas] ./ (scc_gas_pulse_size_conversions[gas] .* options.pulse_size)
    include_slr = base[:DamageAggregator, :include_slr]

    if include_slr
        ciam_mds = _compute_ciam_marginal_damages(base, marginal, gas, ciam_base, ciam_modified, segment_fingerprints; CIAM_foresight=options.CIAM_foresight, CIAM_GDPcap=options.CIAM_GDPcap, pulse_size=options.pulse_size) # NamedTuple with globe and domestic
        # zero out the CIAM marginal damages from start year (2020) through emissions
        # year - they will be non-zero due to foresight but saved marginal damages
        # should be zeroed out pre-emissions year
        ciam_mds.globe[1:year_index] .= 0.
        ciam_mds.domestic[1:year_index] .= 0.
    end

    main_mds = (damages_marginal .- damages_base) .* gas_units_multiplier
    slr_mds = include_slr ? ciam_mds.globe : fill(0., length(_model_years))
    total_mds = main_mds .+ slr_mds

    if options.compute_domestic_values
        main_mds_domestic = (damages_marginal_domestic .- damages_base_domestic) .* gas_units_multiplier
        slr_mds_domestic = include_slr ? ciam_mds.domestic : fill(0., length(_model_years))
        total_mds_domestic = main_mds_domestic .+ slr_mds_domestic
    end

    if options.compute_sectoral_values
        cromar_mortality_mds    = (marginal[:DamageAggregator, :cromar_mortality_damage] .- base[:DamageAggregator, :cromar_mortality_damage]) .* gas_units_multiplier
        agriculture_mds         = (marginal[:DamageAggregator, :agriculture_damage] .- base[:DamageAggregator, :agriculture_damage]) .* gas_units_multiplier
        energy_mds              = (marginal[:DamageAggregator, :energy_damage] .- base[:DamageAggregator, :energy_damage]) .* gas_units_multiplier 
    
        if options.compute_domestic_values
            cromar_mortality_mds_domestic    = (marginal[:DamageAggregator, :cromar_mortality_damage_domestic] .- base[:DamageAggregator, :cromar_mortality_damage_domestic]) .* gas_units_multiplier
            agriculture_mds_domestic         = (marginal[:DamageAggregator, :agriculture_damage_domestic] .- base[:DamageAggregator, :agriculture_damage_domestic]) .* gas_units_multiplier
            energy_mds_domestic              = (marginal[:DamageAggregator, :energy_damage_domestic] .- base[:DamageAggregator, :energy_damage_domestic]) .* gas_units_multiplier 
        end
    end

    # Save marginal damages
    if options.save_md

        # global
        md_values[(region=:globe, sector=:total)][trialnum, :] = total_mds[_damages_idxs]
        if options.compute_sectoral_values
            md_values[(region=:globe, sector=:cromar_mortality)][trialnum, :]   = cromar_mortality_mds[_damages_idxs]
            md_values[(region=:globe, sector=:agriculture)][trialnum, :]        = agriculture_mds[_damages_idxs]
            md_values[(region=:globe, sector=:energy)][trialnum, :]             = energy_mds[_damages_idxs]
            md_values[(region=:globe, sector=:slr)][trialnum, :]                = slr_mds[_damages_idxs]
        end

        # domestic
        if options.compute_domestic_values
            md_values[(region=:domestic, sector=:total)][trialnum, :] = total_mds_domestic[_damages_idxs]
            if options.compute_sectoral_values
                md_values[(region=:domestic, sector=:cromar_mortality)][trialnum, :]   = cromar_mortality_mds_domestic[_damages_idxs]
                md_values[(region=:domestic, sector=:agriculture)][trialnum, :]        = agriculture_mds_domestic[_damages_idxs]
                md_values[(region=:domestic, sector=:energy)][trialnum, :]             = energy_mds_domestic[_damages_idxs]
                md_values[(region=:domestic, sector=:slr)][trialnum, :]                = slr_mds_domestic[_damages_idxs]
            end
        end
    end

    # Save slr damages
    if options.save_slr_damages

        # get a dummy ciam model to be sure to accurately assign segment names to 
        # segment level damages
        m = MimiGIVE.get_model()
        m_ciam, ~ = MimiGIVE.get_ciam(m)

        if include_slr

            # global
            slr_damages[:base][trialnum,:] = ciam_mds.damages_base[_damages_idxs]
            slr_damages[:modified][trialnum,:] = ciam_mds.damages_modified[_damages_idxs]
            slr_damages[:base_lim_cnt][trialnum,:,:] = ciam_mds.base_lim_cnt
            slr_damages[:modified_lim_cnt][trialnum,:,:] = ciam_mds.modified_lim_cnt
            slr_damages[:base_segments_2100][trialnum, :] = ciam_mds.damages_base_segments_2100

            # domestic - these Dictionary entries will only exist if we are computing
            # domestic values
            if options.compute_domestic_values
                slr_damages[:base_domestic][trialnum,:] = ciam_mds.damages_base_domestic[_damages_idxs]
                slr_damages[:modified_domestic][trialnum,:] = ciam_mds.damages_modified_domestic[_damages_idxs]
            end

        else

            # global
            slr_damages[:base][trialnum,:] .= 0.
            slr_damages[:modified][trialnum,:] .= 0.
            slr_damages[:base_lim_cnt][trialnum,:,:] .= 0.
            slr_damages[:modified_lim_cnt][trialnum,:,:] .= 0.
            slr_damages[:base_segments_2100][trialnum, :] .= 0.

            # domestic - these Dictionary entries will only exist if we are computing
            # domestic values
            if options.compute_domestic_values
                slr_damages[:base_domestic][trialnum,:] .= 0.
                slr_damages[:modified_domestic][trialnum,:] .= 0.
            end
        end
    end

    # Get per capita consumption
    # We don't care about units here because we are only going to use ratios
    cpc = base[:global_netconsumption, :net_cpc]
    
    # Save per capita consumption
    if options.save_cpc
        cpc_values[(region=:globe, sector=:total)][trialnum, :] = cpc[_damages_idxs]
    end

    # Calculate the SCC for each discount rate
    for dr in discount_rates
        df = [((cpc[year_index]/cpc[i])^dr.eta * 1/(1+dr.prtp)^(t-year) for (i,t) in enumerate(_model_years) if year<=t<=last_year)...]
        if options.certainty_equivalent
            df_ce = [((1. / cpc[i])^dr.eta * 1/(1+dr.prtp)^(t-year) for (i,t) in enumerate(_model_years) if year<=t<=last_year)...] # only used if optionas.certainty_equivalent=true
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

            if options.certainty_equivalent
                intermediate_ce_scc = sum(df_ce .* cromar_mortality_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
    
                intermediate_ce_scc = sum(df_ce .* agriculture_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
    
                intermediate_ce_scc = sum(df_ce .* energy_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
    
                intermediate_ce_scc = sum(df_ce .* slr_mds[year_index:last_year_index])
                intermediate_ce_scc_values[(region=:globe, sector=:slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc    
            end

            # sectoral domestic (region=:domestic)
            if options.compute_domestic_values

                scc = sum(df .* cromar_mortality_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector= :cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc
    
                scc = sum(df .* agriculture_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector= :agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc
    
                scc = sum(df .* energy_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector= :energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc
    
                scc = sum(df .* slr_mds_domestic[year_index:last_year_index])
                scc_values[(region=:domestic, sector= :slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = scc
    
                if options.certainty_equivalent
                    intermediate_ce_scc = sum(df_ce .* cromar_mortality_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector= :cromar_mortality, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
        
                    intermediate_ce_scc = sum(df_ce .* agriculture_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector= :agriculture, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
        
                    intermediate_ce_scc = sum(df_ce .* energy_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector= :energy, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc
        
                    intermediate_ce_scc = sum(df_ce .* slr_mds_domestic[year_index:last_year_index])
                    intermediate_ce_scc_values[(region=:domestic, sector= :slr, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta)][trialnum] = intermediate_ce_scc    
                end
            end
        end
    end
end

function _compute_scc_mcs(mm::MarginalModel, 
                            n; 
                            year::Int, 
                            last_year::Int, 
                            discount_rates, 
                            certainty_equivalent::Bool,
                            fair_parameter_set::Symbol = :random,
                            fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
                            rffsp_sampling::Symbol = :random,
                            rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
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

    socioeconomics_module = _get_module_name(mm.base, :Socioeconomic)
    if socioeconomics_module == :MimiSSPs
        socioeconomics_source = :SSP
    elseif socioeconomics_module == :MimiRFFSPs
        socioeconomics_source = :RFF
    end

    mcs = get_mcs(n; 
                    socioeconomics_source=socioeconomics_source, 
                    mcs_years = _model_years, 
                    fair_parameter_set = fair_parameter_set,
                    fair_parameter_set_ids = fair_parameter_set_ids,
                    rffsp_sampling = rffsp_sampling,
                    rffsp_sampling_ids = rffsp_sampling_ids,
                    save_list = save_list
                )
    
    if post_mcs_creation_function!==nothing
        post_mcs_creation_function(mcs)
    end

    regions = compute_domestic_values ? [:globe, :domestic] : [:globe]
    sectors = compute_sectoral_values ? [:total,  :cromar_mortality, :agriculture, :energy, :slr] : [:total]

    scc_values = Dict((region=r, sector=s, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta) => Vector{Float64}(undef, n) for dr in discount_rates, r in regions, s in sectors)
    intermediate_ce_scc_values = certainty_equivalent ? Dict((region=r, sector=s, dr_label=dr.label, prtp=dr.prtp, eta=dr.eta) => Vector{Float64}(undef, n) for dr in discount_rates, r in regions, s in sectors) : nothing
    md_values = save_md ? Dict((region=r, sector=s) => Array{Float64}(undef, n, length(_damages_years)) for r in regions, s in sectors) : nothing
    cpc_values = save_cpc ? Dict((region=r, sector=s) => Array{Float64}(undef, n, length(_damages_years)) for r in [:globe], s in [:total]) : nothing # just global and total for now
    
    if save_slr_damages

        # global
        slr_damages = Dict(
            :base               => Array{Float64}(undef, n, length(_damages_years)),
            :modified           => Array{Float64}(undef, n, length(_damages_years)),
            :base_lim_cnt       => Array{Float64}(undef, n, length(_damages_years), 145), # 145 CIAM countries
            :modified_lim_cnt   => Array{Float64}(undef, n, length(_damages_years), 145), # 145 CIAM countries
            :base_segments_2100 => Array{Float64}(undef, n, 11835) # 11,835 segments
        )

        # domestic
        # optionally add arrays to hold the domestic base and modified damages
        if compute_domestic_values
            slr_damages[:base_domestic] = Array{Float64}(undef, n, length(_damages_years))
            slr_damages[:modified_domestic] = Array{Float64}(undef, n, length(_damages_years))
        end
    else
        slr_damages = nothing
    end

    ciam_base, segment_fingerprints = get_ciam(mm.base)
    ciam_modified, _ = get_ciam(mm.base)

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
                        post_trial_func = post_trial_func,
                        results_in_memory = false,
                        results_output_dir = "$output_dir/results"
                    )

    # unpack the payload object
    scc_values, intermediate_ce_scc_values, md_values, cpc_values, slr_damages, year, last_year, discount_rates, gas, ciam_base, ciam_modified, segment_fingerprints, options = Mimi.payload2(sim_results)
    
    # Write out the slr damages to disk in the same place that variables from the save_list would be written out
    if save_slr_damages
        isdir("$output_dir/results/model_1") || mkpath("$output_dir/results/model_1")
        isdir("$output_dir/results/model_2") || mkpath("$output_dir/results/model_2")

        # global 
        df = DataFrame(slr_damages[:base], :auto) |> 
            i -> rename!(i, Symbol.(_damages_years)) |> 
            i -> insertcols!(i, 1, :trial => 1:n) |> 
            i -> stack(i, Not(:trial)) |>
            i -> rename!(i, [:trial, :time, :slr_damages]) |>
            save("$output_dir/results/model_1/slr_damages.csv")

        df = DataFrame(slr_damages[:modified], :auto) |> 
            i -> rename!(i, Symbol.(_damages_years)) |> 
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
                    i -> rename!(i, Symbol.(_damages_years)) |> 
                    i -> insertcols!(i, 1, :trial => 1:n) |> 
                    i -> stack(i, Not(:trial)) |>
                    i -> rename!(i, [:trial, :time, :slr_damages_domestic]) |>
                    save("$output_dir/results/model_1/slr_damages_domestic.csv")

                df = DataFrame(slr_damages[:modified_domestic], :auto) |> 
                    i -> rename!(i, Symbol.(_damages_years)) |> 
                    i -> insertcols!(i, 1, :trial => 1:n) |> 
                    i -> stack(i, Not(:trial)) |>
                    i -> rename!(i, [:trial, :time, :slr_damages_domestic]) |>
                    save("$output_dir/results/model_2/slr_damages_domestic.csv")
        end

        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))

        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))
        df = DataFrame(:trial => [], :time => [], :country => [], :capped_flag => [])
        for trial in 1:n # loop over trials
            trial_df = DataFrame(slr_damages[:base_lim_cnt][trial,:,:], :auto) |>
                i -> rename!(i, ciam_country_names) |>
                i -> insertcols!(i, 1, :time => _damages_years) |> 
                i -> stack(i, Not(:time)) |>
                i -> insertcols!(i, 1, :trial => fill(trial, length(_damages_years) * 145)) |>
                i -> rename!(i, [:trial, :time, :country, :capped_flag]) |>
                i -> @filter(i, _.capped_flag == 1) |>
                DataFrame
            append!(df, trial_df)
        end
        df |> save("$output_dir/results/slr_damages_base_lim_counts.csv")

        df = DataFrame(:trial => [], :time => [], :country => [], :capped_flag => [])
        ciam_country_names = Symbol.(dim_keys(ciam_base, :ciam_country))
        for trial in 1:n # loop over trials
            trial_df = DataFrame(slr_damages[:modified_lim_cnt][trial,:,:], :auto) |>
                i -> rename!(i, ciam_country_names) |>
                i -> insertcols!(i, 1, :time => _damages_years) |> 
                i -> stack(i, Not(:time)) |>
                i -> insertcols!(i, 1, :trial => fill(trial, length(_damages_years) * 145)) |>
                i -> rename!(i, [:trial, :time, :country, :capped_flag]) |>
                i -> @filter(i, _.capped_flag == 1) |>
                DataFrame
            append!(df, trial_df)
        end
        df |> save("$output_dir/results/slr_damages_modified_lim_counts.csv")
    end

    expected_mu_in_year_of_emission = Dict()

    if certainty_equivalent
        year_index = findfirst(isequal(year), _damages_years)
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
    for (k,v) in scc_values
        if certainty_equivalent
            result[:scc][k] = (
                expected_scc = mean(v),
                se_expected_scc = std(v) / sqrt(n),
                ce_scc = mean(intermediate_ce_scc_values[k]) ./ expected_mu_in_year_of_emission[k],
                ce_sccs= intermediate_ce_scc_values[k] ./ expected_mu_in_year_of_emission[k],
                sccs = v,                
            )
        else
            result[:scc][k] = (
                expected_scc = mean(v),
                se_expected_scc = std(v) / sqrt(n),
                sccs = v
            )
        end
    end

    # add a :mds dictionary, where key value pairs (k,v) are NamedTuples with keys(region, sector) => values are (n x 281 (2020:2300)) matrices
    if save_md
        result[:mds] = Dict()
        for (k,v) in md_values
            result[:mds][k] = v
        end
    end

    # add a :cpc dictionary, where key value pairs (k,v) are NamedTuples with keys(region, sector) => values are (n x 281 (2020:2300)) matrices
    if save_cpc
        result[:cpc] = Dict()
        for (k,v) in cpc_values
            result[:cpc][k] = v
        end
    end

    return result
end

function _compute_ciam_marginal_damages(base, modified, gas, ciam_base, ciam_modified, segment_fingerprints; CIAM_foresight, CIAM_GDPcap, pulse_size)
    update_ciam!(ciam_base, base, segment_fingerprints)
    update_ciam!(ciam_modified, modified, segment_fingerprints)

    run(ciam_base)
    run(ciam_modified)

    # Adjust to use perfect foresight if CIAM_foresight == :perfect
    if CIAM_foresight == :perfect
        OptimalCost_base = compute_PerfectForesight_OptimalCosts(ciam_base)
        OptimalCost_modified = compute_PerfectForesight_OptimalCosts(ciam_modified)
    elseif CIAM_foresight == :limited
        OptimalCost_base = ciam_base[:slrcost, :OptimalCost]
        OptimalCost_modified = ciam_modified[:slrcost, :OptimalCost]
    else
        error("CIAM_foresight must be either :limited or :perfect.")
    end

    # Aggregate to Country-Level Damages
    
    # Obtain a key mapping segment ids to ciam country ids, both of which
    # line up with the orders of dim_keys of ciam_base
    xsc = ciam_base[:slrcost, :xsc]::Dict{Int, Tuple{Int, Int, Int}} 
    ciam_country_mapping = DataFrame(:segment_id => collect(keys(xsc)), :ciam_country_id => first.(collect(values(xsc))))
    
    num_ciam_countries = length(dim_keys(ciam_base, :ciam_country))

    OptimalCost_base_country = Array{Float64}(undef, length(_damages_years), num_ciam_countries)
    OptimalCost_modified_country = Array{Float64}(undef, length(_damages_years), num_ciam_countries)

    for country in 1:num_ciam_countries # 145 consecutive Region IDs mapping to the 145 countries in ciam_base dimension ciam_country

        rows = [findall(i -> i == country, ciam_country_mapping.ciam_country_id)...] # rows of the mapping DataFrame that have this ciam country
        matching_segment_ids = [ciam_country_mapping.segment_id[rows]...] # the actual segment IDs that map to this ciam country

        base_damages = sum(view(OptimalCost_base, :, matching_segment_ids), dims=2)
        OptimalCost_base_country[:, country] = [repeat(base_damages[1:end-1], inner=10); base_damages[end]] # repeat to annual from decadal

        modified_damages = sum(view(OptimalCost_modified, :, matching_segment_ids), dims=2)
        OptimalCost_modified_country[:, country] = [repeat(modified_damages[1:end-1], inner=10); modified_damages[end]] # repeat to annual from decadal
    end

    # Limit Country-Level Sea Level Rise Damages to Country-Level GDP

    if CIAM_GDPcap
        # Obtain annual country-level GDP, select 2020:2300 and CIAM countries, convert from $2005 to $2010 to match CIAM
        gdp = base[:Socioeconomic, :gdp][_damages_idxs, indexin(dim_keys(ciam_base, :ciam_country), dim_keys(base, :country))] .* 1 / pricelevel_2010_to_2005

        base_lim_cnt = Int64.(OptimalCost_base_country .> gdp)
        modified_lim_cnt = Int64.(OptimalCost_modified_country .> gdp)

        OptimalCost_base_country = min.(OptimalCost_base_country, gdp)
        OptimalCost_modified_country = min.(OptimalCost_modified_country, gdp)
    else
        base_lim_cnt = fill(0., length(_damages_years), num_ciam_countries)
        modified_lim_cnt = fill(0., length(_damages_years), num_ciam_countries)
    end

    # domestic
    damages_base_domestic = vec(sum(OptimalCost_base_country[:,138],dims=2)) .* pricelevel_2010_to_2005 # Unit of CIAM is billion USD $2010, convert to billion USD $2005
    damages_modified_domestic = vec(sum(OptimalCost_modified_country[:,138],dims=2)) .* pricelevel_2010_to_2005 # Unit of CIAM is billion USD $2010, convert to billion USD $2005

    damages_marginal_domestic = (damages_modified_domestic .- damages_base_domestic) .* scc_gas_molecular_conversions[gas] ./ (scc_gas_pulse_size_conversions[gas] .* pulse_size) # adjust for the (1) molecular mass and (2) pulse size
    damages_marginal_domestic = damages_marginal_domestic .* 1e9  # Unit at this point is billion USD $2005, we convert to just USD here
    
    # global
    damages_base = vec(sum(OptimalCost_base_country,dims=2)) .* pricelevel_2010_to_2005 # Unit of CIAM is billion USD $2010, convert to billion USD $2005
    damages_modified = vec(sum(OptimalCost_modified_country,dims=2)) .* pricelevel_2010_to_2005 # Unit of CIAM is billion USD $2010, convert to billion USD $2005

    damages_marginal = (damages_modified .- damages_base) .* scc_gas_molecular_conversions[gas] ./ (scc_gas_pulse_size_conversions[gas] .* pulse_size) # adjust for the (1) molecular mass and (2) pulse size
    damages_marginal = damages_marginal .* 1e9 # Unit at this point is billion USD $2005, we convert to just USD here

    # CIAM starts in 2020 so pad with zeros at the beginning
    return (globe               = [fill(0., 2020 - _model_years[1]); damages_marginal], # billion USD $2005
            domestic            = [fill(0., 2020 - _model_years[1]); damages_marginal_domestic], # billion USD $2005
            damages_base        = [fill(0., 2020 - _model_years[1]); damages_base], # billion USD $2005
            damages_modified    = [fill(0., 2020 - _model_years[1]); damages_modified], # billion USD $2005
            damages_base_domestic       = [fill(0., 2020 - _model_years[1]); damages_base_domestic], # billion USD $2005
            damages_modified_domestic = [fill(0., 2020 - _model_years[1]); damages_modified_domestic], # billion USD $2005
            base_lim_cnt        = base_lim_cnt, # 2020:2300 x countries
            modified_lim_cnt    = modified_lim_cnt, # 2020:2300 x countries
            damages_base_segments_2100   = OptimalCost_base[9, :] .* pricelevel_2010_to_2005 # billion USD $2005, 2100 is index 9 in 2020:10:2300, this is uncapped segment-level baseline damages in 2100
    )
end

"""
    get_marginal_model(m::Model; year::Union{Int, Nothing} = nothing, gas::Symbol, pulse_size::Float64)

Creates a Mimi MarginalModel where the provided m is the base model, and the marginal model has additional emissions in year `year`.
The marginal model will have an additional `pulse_size` of emissions in the specified `year` for gas `gas`, which will be in 
units of GtC for CO2, MtN2 for N2O, and MtCH4 for CH4. If no Model m is provided, the default model from MimiGIVE.get_model() 
is used as the base model.
"""
function get_marginal_model(m::Model; year::Union{Int, Nothing} = nothing, gas::Symbol, pulse_size::Float64)
    year === nothing ? error("Must specify an emission year. Try `get_marginal_model(m, year=2020)`.") : nothing
    !(year in _model_years) ? error("Cannot add marginal emissions in $year, year must be within the model's time index $_model_years.") : nothing

    # note here that the pulse size will be used as the `delta` parameter for 
    # the `MarginalModel` and thus allow computation of the SCC to return units of
    # dollars per ton, as long as `pulse_size` is interpreted as baseline units
    # of the given gas, which is units of GtC for CO2, MtN2 for N2O, and MtCH4 for CH4.
    mm = create_marginal_model(m, scc_gas_pulse_size_conversions[gas])
    add_marginal_emissions!(mm.modified, year, gas, pulse_size)

    return mm

end

"""
    add_marginal_emissions!(m::Model, year::Int, gas::Symbol, pulse_size::Float64)

Adds a marginal emission component to year m which adds the pulse_size of additional emissions in the specified `year` for gas `gas`, 
which will be in units of GtC for CO2, MtN2 for N2O, and MtCH4 for CH4.
"""
function add_marginal_emissions!(m::Model, year::Int, gas::Symbol, pulse_size::Float64) 

    time = Mimi.dim_keys(m, :time)
    pulse_year_index = findfirst(i -> i == year, time)

    hfc_list = [:HFC23, :HFC32, :HFC43_10, :HFC125, :HFC134a, :HFC143a, :HFC227ea, :HFC245fa]

    if gas == :CO2

        add_comp!(m, Mimi.adder, :marginalemission, before=:co2_cycle)

        addem = zeros(length(time))
        addem[pulse_year_index] = pulse_size     # GtC in this year

        set_param!(m, :marginalemission, :add, addem)

        connect_param!(m, :marginalemission => :input, :co2_emissions_identity => :output_co2)
        connect_param!(m, :co2_cycle => :E_co2, :marginalemission => :output)

    elseif gas == :CH4

        add_comp!(m, Mimi.adder, :marginalemission, before=:ch4_cycle)

        addem = zeros(length(time))
        addem[pulse_year_index] = pulse_size     # MtCH4 in this year

        set_param!(m, :marginalemission, :add, addem)

        connect_param!(m, :marginalemission => :input, :ch4_emissions_identity => :output_ch4)
        connect_param!(m, :ch4_cycle => :fossil_emiss_CH₄, :marginalemission => :output)

    elseif gas == :N2O
        
        add_comp!(m, Mimi.adder, :marginalemission, before=:n2o_cycle)

        addem = zeros(length(time))
        addem[pulse_year_index] = pulse_size     # MtN2 in this year

        set_param!(m, :marginalemission, :add, addem)

        connect_param!(m, :marginalemission => :input, :n2o_emissions_identity => :output_n2o)
        connect_param!(m, :n2o_cycle => :fossil_emiss_N₂O, :marginalemission => :output)

    elseif gas in hfc_list

        # get gas index
        other_ghg_gases = Mimi.dim_keys(m, :other_ghg)
        gas_index = findfirst(i -> i == gas, Symbol.(other_ghg_gases))

        # perturb hfc emissions

        # For now this will return :emiss_other_ghg because it is treated as a 
        # shared parameter in MimiFAIRv1_6_2, and thus also in this model, but this
        # line keeps us robust if it becomes an unshared parameter.
        model_param_name = Mimi.get_model_param_name(m, :other_ghg_cycles, :emiss_other_ghg)

        # Obtain the base emissions values from the model - the following line 
        # allows us to do so without running the model. If we had run the model
        # we can use deepcopy(m[:other_ghg_cycles, :emiss_other_ghg])
        new_emissions = deepcopy(Mimi.model_param(m, model_param_name).values.data)

        # update emissions parameter with a pulse
        new_emissions[pulse_year_index, gas_index] +=  1.0 # add 1 kt hfc
        update_param!(m, :emiss_other_ghg, new_emissions)
        
    else
        error("Gas `" + gas + "` is not supported.")
    end
end
