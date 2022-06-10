using Mimi, CSVFiles, DataFrames, Query, StatsBase, XLSX, Interpolations, DelimitedFiles, Distributions

function get_ciam(m_give::Mimi.Model)
    
    # --------------------------------------------------------------------------
    # Model Parameters

    # load the list of countries to use for CIAM, these are CIAM countries intersected
    # with available CIAM components
    ciam_countries  = (load(joinpath(@__DIR__, "..", "data", "Dimension_ciam_countries.csv")) |> @select(:CountryISO) |> DataFrame |> Matrix)[:]
    countries = Mimi.dim_keys(m_give, :country)
    segments = Mimi.dim_keys(m_give, :segments)

    # get segment fingerprint information
    segment_fingerprints = load(joinpath(@__DIR__, "../data/CIAM/segment_fingerprints.csv"))  |>
        DataFrame |>
        @filter(_.rgn in ciam_countries) |> # only the segments in the coastal countries we are using
        DataFrame
    segments != segment_fingerprints.segments && error("The segments in segment_fingerprints key need to match the segments in m_give.")

    # time parameters
    tstep = 10 # this is assumed within the slrcost component -- DO NOT CHANGE
    period_length = 50
    start_year = 2020
    end_year = 2300
    times = start_year:tstep:end_year
    adaptPers = convert.(Int64, indexin(unique([start_year, start_year:period_length:end_year..., end_year]), times))

    # --------------------------------------------------------------------------
    # Model Construction

    m = Model()
    set_dimension!(m, :time, length(times))
    set_dimension!(m, :ciam_country, ciam_countries) # countries used in CIAM that overlap with socioeconomics 
    set_dimension!(m, :segments, segments) # segments in the ciam_countries
    set_dimension!(m, :adaptPers, length(adaptPers))
    
    # Add CIAM components
    add_comp!(m, MimiCIAM.slrcost)

    # --------------------------------------------------------------------------
    # The Rest of the Parameters
    
    rgns, segs, ciam_params = get_ciam_params(;first = start_year, tstep = tstep, last = end_year, 
                    adaptation_firsts = adaptPers, ciam_countries = ciam_countries,
                    xsc_params_path = joinpath(@__DIR__,"..","data","CIAM", "xsc_ciam_countries.csv")
                    )
    
	# Check Dimensions
    Mimi.dim_keys(m, :ciam_country) != rgns && error("The countries in xsc key need to match the segments in m_give.")
    Mimi.dim_keys(m, :segments) != segs && error("The segments in xsc key need to match the segments in m_give.")

    for (k,v) in ciam_params
        # these are parameters we don't need to set, the correct one for the run
        # is held in "surgeexposure"
        if !(k in ["surgeexposure_dc-gtsr", "surgeexposure_gtsr"]) 
            update_param!(m, :slrcost, Symbol(k), v) 
        end
    end

    # Set dummy variables

    update_param!(m, :slrcost, :lslr, zeros(length(times), length(segments))) # local sea level rise in meters
    update_param!(m, :slrcost, :pop, zeros(length(times), length(ciam_countries))) # population in millions
    update_param!(m, :slrcost, :ypcc, zeros(length(times), length(ciam_countries))) # ypcc in USD in $2010/yr/person
    update_param!(m, :slrcost, :vsl_ciam_country, zeros(length(times), length(ciam_countries)))

    return m, segment_fingerprints
end
 
function update_ciam!(m, m_give, segment_fingerprints)
    # time parameters
    tstep = 10 # this is assumed within the slrcost component -- DO NOT CHANGE
    start_year = 2020
    end_year = 2300
    normalization_year = 2000 # normalize sea level rise to 2000
    times = start_year:tstep:end_year

    # --------------------------------------------------------------------------
    # Parameters from GIVE Model m_give

    # Indices from GIVE Model m_give to slrcost
    m_give_years = Mimi.dim_keys(m_give, :time)
    time_idxs = convert.(Int64, indexin(start_year:tstep:end_year, m_give_years))
    idx_2000 = findfirst(i -> i == normalization_year, m_give_years)
    country_idxs = convert.(Int64, indexin(Mimi.dim_keys(m, :ciam_country), Mimi.dim_keys(m_give, :country)))

    segments = Mimi.dim_keys(m, :segments)

	# Downscale GMSL to LSL for CIAM segemnts
    lslr = zeros(length(times), length(segments))
    for i in 1:length(times)
        lslr[i,:] = segment_fingerprints.fpGSIC_loc .*  (m_give[:glaciers_small_icecaps, :gsic_sea_level][time_idxs][i] .-  m_give[:glaciers_small_icecaps, :gsic_sea_level][idx_2000]) +
                    segment_fingerprints.fpGIS_loc  .*  (m_give[:greenland_icesheet, :greenland_sea_level][time_idxs][i] .-  m_give[:greenland_icesheet, :greenland_sea_level][idx_2000])+
                    segment_fingerprints.fpAIS_loc  .*  (m_give[:antarctic_icesheet, :ais_sea_level][time_idxs][i] .-  m_give[:antarctic_icesheet, :ais_sea_level][idx_2000])+
                    segment_fingerprints.fpTE_loc   .*  (m_give[:thermal_expansion, :te_sea_level][time_idxs][i] .-  m_give[:thermal_expansion, :te_sea_level][idx_2000])+
                    segment_fingerprints.fpLWS_loc  .*  (m_give[:landwater_storage, :lws_sea_level][time_idxs][i] .-  m_give[:landwater_storage, :lws_sea_level][idx_2000])
    end

    update_param!(m, :slrcost, :lslr, lslr) # local sea level rise in meters

    # Socioeconomics
    #   (1) select the ciam countries from the full set of countries 
    #   (2) convert GDP and Population into Per Capita GDP (2010\$ USD per year) 
    #   per capita) for the slrcost component
    #   (3) get the VSL for each CIAM country
    
    population  = m_give[:Socioeconomic, :population][time_idxs,country_idxs] # millions 
    gdp         = m_give[:Socioeconomic, :gdp][time_idxs,country_idxs] .* 1/pricelevel_2010_to_2005 # billion US $2005/yr -> billion US $2010/yr
    ypcc        = gdp ./ population .* 1000 # USD $2010/yr/person

    update_param!(m, :slrcost, :pop, population) # population in millions
    update_param!(m, :slrcost, :ypcc, ypcc) # ypcc in USD in $2010/yr/person

    # Calculate the VSL
    α       = m_give[:VSL, :α]   # VSL scaling parameter.
    ϵ       = m_give[:VSL, :ϵ]   # Income elasticity of the value of a statistical life.
    y₀      = m_give[:VSL, :y₀]  # Normalization constant.
    # mirror the VSL component which follows the FUND mortality equation:  v.vsl[t,c] = p.α * (p.pc_gdp[t,c] / p.y₀) ^ p.ϵ
    # component expects vsl in millions of US $2010 dollars
    vsl_ciam_country = (α * (ypcc ./ y₀) .^  ϵ) ./ 1e6
    update_param!(m, :slrcost, :vsl_ciam_country, vsl_ciam_country)
end

function compute_PerfectForesight_OptimalCosts_typestable(protect_cost, retreat_cost, no_adapt_cost, ntsteps, nsegments)
    # These are the decision options, each is a permutation of choice and level,
    # that we will allow. Note we ignore ProtectCost0 and RetreatCost0, with idxs 
    # 4 and 5 respectively, because we use allowMaintain = false.
    decision_options = [(label = :ProtectCost10, choice = :ProtectCost, level = 10, idx = 1),
                        (label = :ProtectCost100, choice = :ProtectCost, level = 100, idx = 2),
                        (label = :ProtectCost1000, choice = :ProtectCost, level = 1000, idx = 3),
                        (label = :ProtectCost10000, choice = :ProtectCost, level = 10000, idx = 4),
                        (label = :RetreatCost1, choice = :RetreatCost, level = 1, idx = 1),
                        (label = :RetreatCost10, choice = :RetreatCost, level = 10, idx = 2),
                        (label = :RetreatCost100, choice = :RetreatCost, level = 100, idx = 3),
                        (label = :RetreatCost1000, choice = :RetreatCost, level = 1000, idx = 4),
                        (label = :RetreatCost10000, choice = :RetreatCost, level = 10000, idx = 5),
                        (label = :NoAdaptCost0, choice = :NoAdaptCost, level = 1, idx = 1)
                    ]
    noptions = length(decision_options)

    # this will hold the optimal costs for each segment after considering Perfect Foresight
    optimal_costs = Array{Float64}(undef, ntsteps, nsegments) 

    # Preallocate this array and reuse for each segment
    npv = Vector{Float64}(undef, noptions)

    # Precompute discount factor
    df = [(1.04)^((1-t)*10) for t in 1:ntsteps]

    # loop over segments finding the optimal decision for each in light of perfect foresight
    # NPV and filling in the optimal costs with the undiscounted costs for that decision
    for segment in 1:nsegments
        npv[1:4] .= (sum(protect_cost[t,segment,level] * df[t] * 10 for t in 1:ntsteps) for level in 1:4) # remove the Maintain level (allowMaintain = false for these runs) which is index 5
        npv[5:9] .= (sum(retreat_cost[t,segment,level] * df[t] * 10 for t in 1:ntsteps) for level in 1:5) # remove the Maintain level (allowMaintain = false for these runs) which is index 6
        npv[10]  = sum(no_adapt_cost[t,segment] * df[t] * 10 for t in 1:ntsteps)

        optimal_decision = decision_options[findmin(npv)[2]]
        if optimal_decision.choice == :ProtectCost
            optimal_costs[:, segment] .= view(protect_cost, :, segment, optimal_decision.idx)
        elseif optimal_decision.choice == :RetreatCost
            optimal_costs[:, segment] .= view(retreat_cost, :, segment, optimal_decision.idx)
        elseif optimal_decision.choice == :NoAdaptCost
            optimal_costs[:, segment] .= view(no_adapt_cost, :, segment, optimal_decision.idx)
        else
            error("Unknown option.")
        end
    end
    
    return optimal_costs
end

# NPV foresight correction
# This correction accounts for the fact that the new version of CIAM considers NPV 
# over the current adaptation period (50 years), whereas the previous 
# GAMS version assumes NPV is known across the entire model time horizon (2000-2100, 
# for example).
function compute_PerfectForesight_OptimalCosts(m::Mimi.ModelInstance)
    ntsteps = length(Mimi.dim_keys(m, :time))
    nsegments = length(Mimi.dim_keys(m, :segments))

    optimal_costs = compute_PerfectForesight_OptimalCosts_typestable(
        m[:slrcost, :ProtectCost],
        m[:slrcost, :RetreatCost],
        m[:slrcost, :NoAdaptCost],
        ntsteps,
        nsegments,
    )

    return optimal_costs
end
