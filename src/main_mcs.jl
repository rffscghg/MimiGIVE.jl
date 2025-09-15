using Distributions, Dates, Mimi, CSVFiles, DataFrames, MimiMooreEtAlAgricultureImpacts, StatsBase
import Mimi: SampleStore, add_RV!, add_transform!, add_save!

"""
    get_mcs(trials; 
            socioeconomics_source::Symbol = :RFF, 
            mcs_years = 1750:2300, 
            fair_parameter_set::Symbol = :random,
            fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
            rffsp_sampling::Symbol = :random,
            rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
            save_list::Vector = [],
            Agriculture_gtap::String = "midDF"
        )

Return a Monte Carlo Simulation definition of type Mimi.SimulationDefinition that
holds all random variables and distributions, as assigned to model component/parameter
pairs, that will be used in a Monte Carlo Simulation. 

- `trials` (required) - number of trials to be run, used for presampling
- `socioeconomics_source` (default :RFF) - which source the Socioeconomics component uses
- `fair_parameter_set` (default :random) - :random means FAIR mcs samples will be 
        chosen randomly from the provided sets, while :deterministic means they will 
        be  based on the provided vector of to `fair_parameter_set_ids` keyword argument. 
- `fair_parameter_set_ids` (default nothing) - if `fair_parameter_set` is set 
        to :deterministic, this `n` element vector provides the fair parameter set ids 
        that will be run, otherwise it is set to `nothing` and ignored.
- `rffsp_sampling` (default :random) - which sampling strategy to use for the RFF 
        SPs, :random means RFF SPs will be chosen randomly, while :deterministic means they 
        will be based on the provided vector of to `rffsp_sampling_ids` keyword argument. 
- `rffsp_sampling_ids` (default nothing) - if `rffsp_sampling` is set to :deterministic, 
        this `n` element vector provides the RFF SP ids that will be run, otherwise it is 
        set to `nothing` and ignored.
- `save_list` (default []) - which parameters and varaibles to save for each trial,
        entered as a vector of Tuples (:component_name, :variable_name)
- Agriculture_gtap (default midDF) - specify the `Agriculture_gtap` input parameter as one of 
        `["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]`, indicating which 
        gtap damage function the component should use. 
"""
function get_mcs(trials;
    socioeconomics_source::Symbol=:RFF,
    mcs_years=1750:2300,
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Vector{Int},Nothing}=nothing,
    rffsp_sampling::Symbol=:random,
    rffsp_sampling_ids::Union{Vector{Int},Nothing}=nothing,
    save_list::Vector=[],
    Agriculture_gtap::String="midDF"
)

    # check some argument conditions
    if fair_parameter_set == :deterministic
        isnothing(fair_parameter_set_ids) && error("If `fair_parameter_set` is :determinsitic, must provide a `fair_parameter_set_ids` vector.")
        length(fair_parameter_set_ids) !== trials && error("The length of the provided `fair_parameter_set_ids` vector must be equal to the number of trials ($trials) run.")
        sum(fair_parameter_set_ids .> 2237) > 0. || sum(fair_parameter_set_ids .< 1) > 0. && error("FAIR parameter set ids must be between 1 and 2237, inclusive.")
    end

    if rffsp_sampling == :deterministic
        isnothing(rffsp_sampling_ids) && error("If `rffsp_sampling` is :determinsitic, must provide a `rffsp_sampling_ids` vector.")
        length(rffsp_sampling_ids) !== trials && error("The length of the provided `rffsp_sampling_ids` vector must be equal to the number of trials ($trials) run.")
        sum(rffsp_sampling_ids .> 10_000) > 0. || sum(rffsp_sampling_ids .< 1) > 0. && error("RFF SP sample ids must be between 1 and 10,000, inclusive.")
    end

    # define the Monte Carlo Simulation and add some simple random variables
    mcs = @defsim begin
        dice2016R2_damage.a2 = Normal(0.00236, 0.00236 / 2) # Nordhaus (2017, PNAS) DICE2016 RV 
    end

    # Howard and Sterner (2017) Damage specification table 2 column 3
    hs_μ_3 = [0.595382733860703; 0.259851128136597]
    hs_σ_3 = [0.0322523508274087 -0.0373892045213768
        -0.0373892045213768 0.063496518648112]
    hs_distribution_3 = MvNormal(hs_μ_3, hs_σ_3)
    hs_coefficients_3 = rand(hs_distribution_3, trials)

    add_RV!(mcs, :rv_hs_damage_t2_base_3, SampleStore(hs_coefficients_3[1, :]))
    add_transform!(mcs, :hs_damage, :t2_base_3, :(=), :rv_hs_damage_t2_base_3)

    add_RV!(mcs, :rv_hs_damage_t2_cat_3, SampleStore(hs_coefficients_3[2, :]))
    add_transform!(mcs, :hs_damage, :t2_cat_3, :(=), :rv_hs_damage_t2_cat_3)

    # Howard and Sterner (2017) Damage specification table 2 column 4
    hs_μ_4 = [0.595382733860703; 0.259851128136597; 0.113324887895228]
    hs_σ_4 = [0.0362838946808348 -0.0420628550865489 0.
        -0.0420628550865489 0.0714335834791260 0.
        0. 0. 0.0157459807497214]
    hs_distribution_4 = MvNormal(hs_μ_4, hs_σ_4)
    hs_coefficients_4 = rand(hs_distribution_4, trials)
    add_RV!(mcs, :rv_hs_damage_t2_base_4, SampleStore(hs_coefficients_4[1, :]))
    add_transform!(mcs, :hs_damage, :t2_base_4, :(=), :rv_hs_damage_t2_base_4)

    add_RV!(mcs, :rv_hs_damage_t2_cat_4, SampleStore(hs_coefficients_4[2, :]))
    add_transform!(mcs, :hs_damage, :t2_cat_4, :(=), :rv_hs_damage_t2_cat_4)

    add_RV!(mcs, :rv_hs_damage_t2_prod_4, SampleStore(hs_coefficients_4[3, :]))
    add_transform!(mcs, :hs_damage, :t2_prod_4, :(=), :rv_hs_damage_t2_prod_4)

    # Howard and Sterner (2017) Damage specification table 2 column 7
    hs_μ_7 = [0.318149737017145; 0.362274271711041]
    hs_σ_7 = [0.00953254601993184 -0.00956576259414058
        -0.00956576259414058 0.00970956896549987]
    hs_distribution_7 = MvNormal(hs_μ_7, hs_σ_7)
    hs_coefficients_7 = rand(hs_distribution_7, trials)

    add_RV!(mcs, :rv_hs_damage_t2_base_7, SampleStore(hs_coefficients_7[1, :]))
    add_transform!(mcs, :hs_damage, :t2_base_7, :(=), :rv_hs_damage_t2_base_7)

    add_RV!(mcs, :rv_hs_damage_t2_cat_7, SampleStore(hs_coefficients_7[2, :]))
    add_transform!(mcs, :hs_damage, :t2_cat_7, :(=), :rv_hs_damage_t2_cat_7)

    # Howard and Sterner (2017) Damage specification table 2 column 8
    hs_μ_8 = [0.318149737017145; 0.362274271711041; 0.398230480262918]
    hs_σ_8 = [0.0104404075456397 -0.0104767876031064 0.
        -0.0104767876031064 0.010634289819357 0.
        0. 0. 0.0563560833680617]
    hs_distribution_8 = MvNormal(hs_μ_8, hs_σ_8)
    hs_coefficients_8 = rand(hs_distribution_8, trials)

    add_RV!(mcs, :rv_hs_damage_t2_base_8, SampleStore(hs_coefficients_8[1, :]))
    add_transform!(mcs, :hs_damage, :t2_base_8, :(=), :rv_hs_damage_t2_base_8)

    add_RV!(mcs, :rv_hs_damage_t2_cat_8, SampleStore(hs_coefficients_8[2, :]))
    add_transform!(mcs, :hs_damage, :t2_cat_8, :(=), :rv_hs_damage_t2_cat_8)

    add_RV!(mcs, :rv_hs_damage_t2_prod_8, SampleStore(hs_coefficients_8[3, :]))
    add_transform!(mcs, :hs_damage, :t2_prod_8, :(=), :rv_hs_damage_t2_prod_8)


    # add the socioeconomics RV if the socioeconomics source is Mimi RFF SPs
    # Use SampleStore for a deterministic RFF SP sampling approach, otherwise
    # use an EmpiricalDistribution across all ids (equal probability is assumed 
    # if probabilities not provided)
    if socioeconomics_source == :RFF
        distrib = rffsp_sampling == :random ? EmpiricalDistribution(collect(1:10_000)) : SampleStore(rffsp_sampling_ids)
        add_RV!(mcs, :socio_id_rv, distrib)
        add_transform!(mcs, :Socioeconomic, :id, :(=), :socio_id_rv)
    end

    #add BRICK random variable - assign one Normally distributed RV per year
    for year in mcs_years
        rv_name = Symbol("rv_landwater_storage_$year")
        add_RV!(mcs, rv_name, Normal(0.0003, 0.00018))
        add_transform!(mcs, :landwater_storage, :lws_random_sample, :(=), rv_name, [year])
    end

    BRICK_parameters = load(joinpath(@__DIR__, "..", "data", "BRICK_posterior_parameters_10k.csv")) |> DataFrame

    BRICK_uncertain_parameters = [
        (source_name=:thermal_s0, comp_name=:thermal_expansion, param_name=:te_s₀),
        (source_name=:thermal_alpha, comp_name=:thermal_expansion, param_name=:te_α), (source_name=:glaciers_v0, comp_name=:glaciers_small_icecaps, param_name=:gsic_v₀),
        (source_name=:glaciers_s0, comp_name=:glaciers_small_icecaps, param_name=:gsic_s₀),
        (source_name=:glaciers_beta0, comp_name=:glaciers_small_icecaps, param_name=:gsic_β₀),
        (source_name=:glaciers_n, comp_name=:glaciers_small_icecaps, param_name=:gsic_n), (source_name=:greenland_v0, comp_name=:greenland_icesheet, param_name=:greenland_v₀),
        (source_name=:greenland_a, comp_name=:greenland_icesheet, param_name=:greenland_a),
        (source_name=:greenland_b, comp_name=:greenland_icesheet, param_name=:greenland_b),
        (source_name=:greenland_alpha, comp_name=:greenland_icesheet, param_name=:greenland_α),
        (source_name=:greenland_beta, comp_name=:greenland_icesheet, param_name=:greenland_β), (source_name=:anto_alpha, comp_name=:antarctic_ocean, param_name=:anto_α),
        (source_name=:anto_beta, comp_name=:antarctic_ocean, param_name=:anto_β), (source_name=:ais_gamma, comp_name=:antarctic_icesheet, param_name=:ais_γ),
        (source_name=:ais_alpha, comp_name=:antarctic_icesheet, param_name=:ais_α),
        (source_name=:ais_mu, comp_name=:antarctic_icesheet, param_name=:ais_μ),
        (source_name=:ais_v, comp_name=:antarctic_icesheet, param_name=:ais_ν),
        (source_name=:ais_precip0, comp_name=:antarctic_icesheet, param_name=:ais_precipitation₀),
        (source_name=:ais_kappa, comp_name=:antarctic_icesheet, param_name=:ais_κ),
        (source_name=:ais_flow0, comp_name=:antarctic_icesheet, param_name=:ais_iceflow₀),
        (source_name=:ais_runoff_height0, comp_name=:antarctic_icesheet, param_name=:ais_runoffline_snowheight₀),
        (source_name=:ais_c, comp_name=:antarctic_icesheet, param_name=:ais_c),
        (source_name=:ais_bedheight0, comp_name=:antarctic_icesheet, param_name=:ais_bedheight₀),
        (source_name=:ais_slope, comp_name=:antarctic_icesheet, param_name=:ais_slope),
        (source_name=:ais_lambda, comp_name=:antarctic_icesheet, param_name=:λ),
        (source_name=:ais_temp_threshold, comp_name=:antarctic_icesheet, param_name=:temperature_threshold),
        (source_name=:antarctic_s0, comp_name=:antarctic_icesheet, param_name=:ais_sea_level₀), # DOUBLE CHECK
    ]

    for p in BRICK_uncertain_parameters
        rv_name = Symbol("rv_brick_$(p.source_name)")
        add_RV!(mcs, rv_name, SampleStore(BRICK_parameters[:, p.source_name]))
        add_transform!(mcs, p.comp_name, p.param_name, :(=), rv_name)
    end

    # add Agriculture mcs over gtap region damage function parameterizations
    ag_sample_stores = MimiMooreEtAlAgricultureImpacts.get_probdists_gtap_df(Agriculture_gtap, trials)

    # If ag sample stores are available for a given Agriculture_gtap damage function 
    # then ag_sample_stores will be a Vector, and otherwise will return a 
    # warning and `nothing`. 
    if !isnothing(ag_sample_stores)
        for coef in [1, 2, 3] # three coefficients defined with an anonymous dimension
            for (i, region) in enumerate(["USA", "CAN", "WEU", "JPK", "ANZ", "EEU", "FSU", "MDE", "CAM", "LAM", "SAS", "SEA", "CHI", "MAF", "SSA", "SIS"]) # fund regions for ag
                rv_name = Symbol("rv_gtap_coef$(coef)_$region")
                add_RV!(mcs, rv_name, ag_sample_stores[i, coef])
                add_transform!(mcs, :Agriculture, :gtap_df, :(=), rv_name, [region, coef])
            end
        end
    end

    # add Cromar uncertainty based on coefficients from Cromar et al.
    cromar_coeffs = load(joinpath(@__DIR__, "..", "data", "CromarMortality_damages_coefficients.csv")) |> DataFrame
    cromar_mapping_raw = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_cromar_mortality_regions.csv")) |> DataFrame

    # Get one random variable per region
    for (i, region) in enumerate(cromar_coeffs[!, "Cromar Region Name"])
        rv_name = Symbol("rv_β_mortality_$(region)")
        add_RV!(mcs, rv_name, Normal(cromar_coeffs[i, "Pooled Beta"], cromar_coeffs[i, "Pooled SE"]))
    end

    # add one transform per country asigning each to the appropriate regional random variable
    for row in 1:size(cromar_mapping_raw, 1)
        rv_name = Symbol("rv_β_mortality_$(cromar_mapping_raw.cromar_region[row])")
        add_transform!(mcs, :CromarMortality, :β_mortality, :(=), rv_name, [cromar_mapping_raw.ISO3[row]])
    end

    # add the FAIR random variables and transforms - note this could be done within
    # the @defsim macro but we use the dictionary to make this less verbose

    # Note that if a parameter component is not included in add_transform!, the 
    # parameters are shared model parameters, and each line will create ONE random 
    # variable and assign all component parameters connected to that shared model 
    # parameter to the value taken on by that random variable

    fair_samples_map, fair_samples = get_fair_mcs_params(trials; fair_parameter_set=fair_parameter_set, fair_parameter_set_ids=fair_parameter_set_ids)
    fair_samples_left = deepcopy(fair_samples) # we know we've added everything when this is empty!

    # add and assign all random variables for single dimensional parameters
    for (k, v) in fair_samples_left
        if size(v, 2) == 1 # one column of values
            rv_name = Symbol("rv_$k")
            add_RV!(mcs, rv_name, SampleStore(fair_samples[k][!, 1]))
            add_transform!(mcs, k, :(=), rv_name)
            delete!(fair_samples_left, k)
        end
    end

    # assign one random variable per year with a unique distribution from fair_samples
    # assume F_solar parameter set defines value starting in 1750 with 361 years total
    for year in 1750:2110
        rv_name = Symbol("rv_F_solar_$year")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:F_solar][!, string(year)]))
        add_transform!(mcs, :F_solar, :(=), rv_name, [year])
    end
    delete!(fair_samples_left, :F_solar)

    # Radiative forcing scaling - one distribution per "other" greenhouse gas, and
    # one per ods

    for gas in names(fair_samples[:scale_other_ghg])
        rv_name = Symbol("rv_scale_other_ghg_$(gas)")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:scale_other_ghg][!, gas]))
        add_transform!(mcs, :scale_other_ghg, :(=), rv_name, [gas])
    end
    delete!(fair_samples_left, :scale_other_ghg)

    for ods in names(fair_samples[:scale_ods])
        rv_name = Symbol("rv_scale_ods_$(ods)")
        add_RV!(mcs, rv_name, SampleStore(fair_samples[:scale_ods][!, ods]))
        add_transform!(mcs, :scale_ods, :(=), rv_name, [ods])
    end
    delete!(fair_samples_left, :scale_ods)

    # ocean_heat_capacity takes an anonymous dim of 2 (deep and mixed, should label 
    # explicilty) - anonymouse dims are named with Ints 1 and 2
    rv_name = Symbol("rv_ocean_heat_capacity_1")
    add_RV!(mcs, rv_name, SampleStore(fair_samples[:ocean_heat_capacity][!, "1"]))
    add_transform!(mcs, :ocean_heat_capacity, :(=), rv_name, [1])

    rv_name = Symbol("rv_ocean_heat_capacity_2")
    add_RV!(mcs, rv_name, SampleStore(fair_samples[:ocean_heat_capacity][!, "2"]))
    add_transform!(mcs, :ocean_heat_capacity, :(=), rv_name, [2])

    delete!(fair_samples_left, :ocean_heat_capacity)

    # check if we've added all FAIR parameters
    isempty(fair_samples_left) ? nothing : error("The following FAIR mcs uncertain parameters has not been added to the simulation: $(keys(fair_samples_left))")

    # add the requested saved variables 
    for i in save_list
        add_save!(mcs, i)
    end

    return mcs
end

"""
    run_mcs(;trials::Int64 = 10000, 
            output_dir::Union{String, Nothing} = nothing, 
            save_trials::Bool = false,
            fair_parameter_set::Symbol = :random,
            fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
            rffsp_sampling::Symbol = :random,
            rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
            m::Mimi.Model = get_model(), 
            save_list::Vector = [],
            results_in_memory::Bool = true,
        )

Return the results of a Monte Carlo Simulation with the defined number of trials
and save data into the `output_dir` folder, optionally also saving trials if 
`save_trials` is set to `true.` If no model is provided, use the default model 
returned by get_model().

- `trials` (default 10,000) - number of trials to be run, used for presampling
- `output_dir` (default constructed folder name) - folder to hold results 
- `save_trials` (default false) - whether to save all random variables for all trials to trials.csv 
- `fair_parameter_set` (default :random) - :random means FAIR mcs samples will be 
        chosen randomly from the provided sets, while :deterministic means they will 
        be  based on the provided vector of to `fair_parameter_set_ids` keyword argument. 
- `fair_parameter_set_ids` - (default nothing) - if `fair_parameter_set` is set 
        to :deterministic, this `n` element vector provides the fair parameter set ids 
        that will be run, otherwise it is set to `nothing` and ignored.
- `rffsp_sampling` (default :random) - which sampling strategy to use for the RFF 
        SPs, :random means RFF SPs will be chosen randomly, while :deterministic means they 
        will be based on the provided vector of to `rffsp_sampling_ids` keyword argument. 
- `rffsp_sampling_ids` - (default nothing) - if `rffsp_sampling` is set to :deterministic, 
        this `n` element vector provides the RFF SP ids that will be run, otherwise it is 
        set to `nothing` and ignored.
- `m` (default get_model()) - the model to run the simulation for
- `save_list` (default []) - which parameters and variables to save for each trial,
        entered as a vector of Tuples (:component_name, :variable_name)
- `results_in_memory` (default true) - this should be turned off if you are running 
        into memory problems, data will be streamed out to disk but not saved in memory 
        to the mcs object
"""
function run_mcs(; trials::Int64=10000,
    output_dir::Union{String,Nothing}=nothing,
    save_trials::Bool=false,
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Vector{Int},Nothing}=nothing,
    rffsp_sampling::Symbol=:random,
    rffsp_sampling_ids::Union{Vector{Int},Nothing}=nothing,
    m::Mimi.Model=get_model(),
    save_list::Vector=[],
    results_in_memory::Bool=true,
)

    m = deepcopy(m) # in the case that an `m` was provided, be careful that we don't modify the original

    trials < 2 && error("Must run `run_mcs` function with a `trials` argument greater than 1 due to a Mimi specification about SampleStores.  TO BE FIXED SOON!")

    # Set up output directories
    output_dir = output_dir === nothing ? joinpath(@__DIR__, "../output/mcs/", "MCS $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$trials") : output_dir
    isdir("$output_dir/results") || mkpath("$output_dir/results")

    trials_output_filename = save_trials ? joinpath("$output_dir/trials.csv") : nothing

    socioeconomics_module = _get_module_name(m, :Socioeconomic)
    if socioeconomics_module == :MimiSSPs
        socioeconomics_source = :SSP
    elseif socioeconomics_module == :MimiRFFSPs
        socioeconomics_source = :RFF
    end

    Agriculture_gtap = _get_mooreag_gtap(m)

    # Get an instance of the mcs
    mcs = get_mcs(trials;
        socioeconomics_source=socioeconomics_source,
        mcs_years=Mimi.time_labels(m),
        fair_parameter_set=fair_parameter_set,
        fair_parameter_set_ids=fair_parameter_set_ids,
        rffsp_sampling=rffsp_sampling,
        rffsp_sampling_ids=rffsp_sampling_ids,
        save_list=save_list,
        Agriculture_gtap=Agriculture_gtap
    )

    # run monte carlo trials
    results = run(mcs,
        m,
        trials;
        trials_output_filename=trials_output_filename,
        results_output_dir="$output_dir/results",
        results_in_memory=results_in_memory
    )

    return results
end

"""
    get_fair_mcs_params(n::Int; 
                        fair_parameter_set::Symbol = :random, 
                        fair_parameter_set_ids::Union{Nothing, Vector{Int}}
                    )
                    
Return the FAIR mcs parameters mapped from parameter name to string name, and a dictionary
using the parameter names as keys and a DataFrame holding the values as a value.
If fair_parameter_set is :random (default), then FAIR mcs samples will be chosen 
randomly from the provided sets. If it set to :deterministic they will be the vector
provided by fair_parameter_set_ids.
"""
function get_fair_mcs_params(n::Int;
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Nothing,Vector{Int}}
)

    # check some argument conditions
    if fair_parameter_set == :deterministic
        isnothing(fair_parameter_set_ids) && error("If `fair_parameter_set` is :determinsitic, must provide a `fair_parameter_set_ids` vector.")
        length(fair_parameter_set_ids) !== n && error("The length of the provided `fair_parameter_set_ids` vector must be equal to the number of trials ($n) run.")
        sum(fair_parameter_set_ids .> 2237) .> 0 || sum(fair_parameter_set_ids .< 1) > 0. && error("FAIR parameter set ids must be between 1 and 2237, inclusive.")
    end

    names_map = get_fair_mcs_params_map()
    params_dict = Dict()

    if fair_parameter_set == :deterministic
        samples = fair_parameter_set_ids
    elseif fair_parameter_set == :random
        samples = sample(collect(1:2237), n, replace=true) # randomly sample n sample sets of 2237 options, with replacement
    end

    for (k, v) in names_map
        values = load(joinpath(@__DIR__, "..", "data", "FAIR_mcs", "fair_mcs_params_$v.csv")) |> DataFrame # load the deterministic set of 2237 parameters
        push!(params_dict, k => values[samples, :])
    end
    return names_map, params_dict

end

"""
    get_fair_mcs_params_map()

Return a dictionary of FAIR elements with the FAIR v1.6.2 parameter name being 
the component, parameter pair and the value being the parameter csv name. 
"""
function get_fair_mcs_params_map()
    return Dict(
        :β_CO => "b_aero_CO",
        :scale_CH₄ => "scale_CH4",
        :F_solar => "F_solar",
        :Ψ_CH₄ => "b_tro3_CH4",
        :scale_N₂O => "scale_N2O",
        :CO₂_pi => "C_pi",
        :deep_ocean_efficacy => "deep_ocean_efficacy",
        :scale_bcsnow => "scale_bcsnow",
        :scale_aerosol_direct_OC => "scale_aerosol_direct_OC",
        :b_SOx => "ghan_params_SOx",
        :feedback => "ozone_feedback",
        :scale_O₃ => "scale_O3",
        :b_POM => "ghan_params_b_POM",
        :r0_co2 => "r0",
        :β_NH3 => "b_aero_NH3",
        :lambda_global => "lambda_global",
        :scale_landuse => "scale_landuse",
        :scale_volcanic => "scale_volcanic",
        :scale_aerosol_direct_SOx => "scale_aerosol_direct_SOx",
        :β_NOx => "b_aero_NOx",
        :Ψ_N₂O => "b_tro3_N2O",
        :ocean_heat_capacity => "ocean_heat_capacity",
        :β_OC => "b_aero_OC",
        :scale_solar => "scale_solar",
        :rC_co2 => "rc",
        :scale_aerosol_direct_BC => "scale_aerosol_direct_BC",
        :scale_CH₄_H₂O => "scale_CH4_H2O",
        :scale_aerosol_indirect => "scale_aerosol_indirect",
        :scale_ods => "scale_ods",
        :Ψ_CO => "b_tro3_CO",
        :scale_aerosol_direct_NOx_NH3 => "scale_aerosol_direct_NOx_NH3",
        :scale_other_ghg => "scale_other_ghg",
        :Ψ_NMVOC => "b_tro3_NMVOC",
        :F2x => "F2x",
        :β_SOx => "b_aero_SOx",
        :β_NMVOC => "b_aero_NMVOC",
        :rT_co2 => "rt",
        :β_BC => "b_aero_BC",
        :scale_CO₂ => "scale_CO2",
        :Ψ_ODS => "b_tro3_ODS",
        :scale_aerosol_direct_CO_NMVOC => "scale_aerosol_direct_CO_NMVOC",
        :Ψ_NOx => "b_tro3_NOx",
        :ocean_heat_exchange => "ocean_heat_exchange",
        :ϕ => "ghan_params_Pi"
    )
end
