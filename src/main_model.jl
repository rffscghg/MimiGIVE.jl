using Mimi, CSVFiles, DataFrames, Query, StatsBase, XLSX, Interpolations, DelimitedFiles, Distributions

"""
    get_model(; Agriculture_gtap::String = "midDF",
                socioeconomics_source::Symbol = :RFF,
                SSP_scenario::Union{Nothing, String} = nothing,       
                RFFSPsample::Union{Nothing, Int} = nothing,
                Agriculture_floor_on_damages::Bool = true,
                Agriculture_ceiling_on_benefits::Bool = false,
                vsl::Symbol= :epa
            )
                
Get a GIVE Model with the given argument Settings

-- Socioeconomic -- 

- socioeconomics_source (default :RFF) - The options are :RFF, which uses data from 
    the RFF socioeconomic projections, or :SSP, which uses data from one of the 
    Shared Socioeconomic Pathways

- SSP_scenario (default to nothing) - This setting is used only if one is using 
    the SSPs as the socioeconomics_source, and the current options are "SSP119", 
    "SSP126", "SSP245", "SSP370", "SSP585", and this will be used as follows.
    See the SSPs component here: https://github.com/anthofflab/MimiSSPs.jl for more information.

    (1) Select the population and GDP trajectories for 2020 through 2300, mapping
        each RCMIP scenario to the SSP (SSP1, 2, 3, 5 respectively)
    
    (2) Choose the ar6 scenario for data from 1750 - 2019 and the RCMIP emissions 
        scenario from the MimiSSPs component to pull Leach et al. RCMIP scenario
        data for 2020 to 2300 for CO2, CH4, and N2O.

    (NOTE) that if the socioeconomics_source is :RFF this will not be consequential 
        and ssp245 will be used for the ar6 data from 1750 - 2019 and trace gases 
        from 2020 onwards, while emissions for CO2, CH4, and N2O will come from
        the MimiRFFSPs component.

- RFFSPsample (default to nothing, which will pull the in MimiRFFSPs) - choose
    the sample for which to run the RFF SP. See the RFFSPs component here: 
    https://github.com/rffscghg/MimiRFFSPs.jl.

-- Agriculture -- 

- Agriculture_gtap (default midDF) - specify the `Agriculture_gtap` input parameter as one of 
    `["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]`, indicating which 
    gtap damage function the component should use. 

 - Agriculture_floor_on_damages (default true) - If `Agriculture_gtap_floor_on_damages` = true, then 
    the agricultural damages (negative values of the `agcost` variable) in each 
    timestep will not be allowed to exceed 100% of the size of the  agricultural 
    sector in each region.

- Agriculture_ceiling_on_benefits (default false) - If `Agriculture_gtap_ceiling_on_benefits` = true, 
    then the agricultural benefits (positive values of the `agcost` variable) in 
    each timestep will not be allowed to exceed 100% of the size of the agricultural 
    sector in each region.

-- Other --

- vsl (default :epa) - specify the soruce of the value of statistical life (VSL) being used in the model

"""
function get_model(; Agriculture_gtap::String="midDF",
    socioeconomics_source::Symbol=:RFF,
    SSP_scenario::Union{Nothing,String}=nothing,
    RFFSPsample::Union{Nothing,Int}=nothing,
    Agriculture_floor_on_damages::Bool=true,
    Agriculture_ceiling_on_benefits::Bool=false,
    vsl::Symbol=:epa
)

    # --------------------------------------------------------------------------
    # MODEL - Check Arguments
    # --------------------------------------------------------------------------    

    if socioeconomics_source == :SSP && isnothing(SSP_scenario)
        error("The socioeconomics_source argument :SSP requires setting a SSP_scenario")
    end

    if socioeconomics_source == :RFF && !isnothing(SSP_scenario)
        @warn("You have set a SSP_scenario to a non-nothing value, but note that setting the socioeconomics_source argument to :RFF means that this will have no effect on the model.")
    end

    # Restrictions on arguments
    socioeconomics_source_options = [:SSP, :RFF]
    socioeconomics_source in socioeconomics_source_options ? nothing : error("The socioeconomics_source must be one of $(socioeconomics_source_options)")
    Agriculture_gtap in MimiMooreEtAlAgricultureImpacts.gtaps ? nothing : error("Unknown GTAP dataframe specification: \"$Agriculture_gtap\". Must be one of the following: $(MimiMooreEtAlAgricultureImpacts.gtaps)")

    SSP_scenario_options = [nothing, "SSP119", "SSP126", "SSP245", "SSP370", "SSP585"]
    SSP_scenario in SSP_scenario_options ? nothing : error("The SSP_scenario must be one of $(SSP_scenario_options)")

    # --------------------------------------------------------------------------
    # MODEL - Model Data and Settings
    # --------------------------------------------------------------------------    

    # dimensions and countries/regions lists
    countries = (load(joinpath(@__DIR__, "..", "data", "Dimension_countries.csv"))|>@select(:CountryISO)|>DataFrame|>Matrix)[:]
    ciam_countries = (load(joinpath(@__DIR__, "..", "data", "Dimension_ciam_countries.csv"))|>@select(:CountryISO)|>DataFrame|>Matrix)[:]
    fund_regions = (load(joinpath(@__DIR__, "..", "data", "Dimension_fund_regions.csv"))|>@select(:fund_region)|>DataFrame|>Matrix)[:]
    gcam_regions = (load(joinpath(@__DIR__, "..", "data", "Dimension_gcam_energy_regions.csv"))|>@select(:gcam_energy_region)|>DataFrame|>Matrix)[:] # not currently a dimension in model
    cromar_regions = (load(joinpath(@__DIR__, "..", "data", "Dimension_cromar_mortality_regions.csv"))|>@select(:cromar_mortality_region)|>DataFrame|>Matrix)[:] # not currently a dimension in model
    domestic_countries = ["USA", "PRI"] # Country ISO3 codes to be accumulated for domestic

    # Create country-region (FUND derived) mapping for Agriculture damage function
    ag_mapping = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_fund_regions.csv")) |> DataFrame
    ag_mapping.ISO3 != countries && error("FUND mapping file ISO3 column must match model countries vector exactly.")
    sort(unique(ag_mapping.fundregion)) != sort(fund_regions) && error("FUND mapping file fund_regions column must match model fund_regions vector exactly (when both are sorted).")
    ag_mapping = ag_mapping.fundregion

    # Create country-region mapping for GCAM energy damage function.
    energy_mapping = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_gcam_energy_regions.csv")) |> DataFrame
    energy_mapping.ISO3 != countries && error("GCAM mapping file ISO3 column must match model countries vector exactly.")
    sort(unique(energy_mapping.gcamregion)) != sort(gcam_regions) && error("GCAM mapping file gcam_regions column must match model gcamregions vector exactly (when both are sorted).")
    energy_mapping = energy_mapping.gcamregion

    # Create country-region mapping for Cromar et al. temperature-mortality damage function.
    cromar_mapping = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_cromar_mortality_regions.csv")) |> DataFrame
    cromar_mapping.ISO3 != countries && error("Cromar mortality mapping file ISO3 column must match model countries vector exactly.")
    sort(unique(cromar_mapping.cromar_region)) != sort(cromar_regions) && error("Cromar mortality mapping file gcam_regions column must match model gcamregions vector exactly (when both are sorted).")
    cromar_mapping = cromar_mapping.cromar_region

    # BRICK Fingerprinting
    segment_fingerprints = load(joinpath(@__DIR__, "../data/CIAM/segment_fingerprints.csv")) |>
                           DataFrame |>
                           @filter(_.rgn in ciam_countries) |> # reduce to the segments in the coastal countries we are using
                           DataFrame

    # get the ar6 forcing scenario to be used for the FAIR model and Mortality component
    if socioeconomics_source == :RFF
        ar6_scenario = "ssp245" # use SSP245 emissions scenario as the basis for trace gases for RFF SP
    elseif socioeconomics_source == :SSP
        ar6_scenario = lowercase(SSP_scenario)
    end

    # Baseline mortality use SSP2 as a proxy for SSP4 and SSP1 as a proxy for 
    # SSP5 per instructions from the literature
    mortality_SSP_map = Dict("SSP1" => "SSP1", "SSP2" => "SSP2", "SSP3" => "SSP3", "SSP4" => "SSP2", "SSP5" => "SSP1")

    # Grab the SSP name from the full scenario ie. SSP2 from SSP245
    if socioeconomics_source == :SSP
        SSP = SSP_scenario[1:4]
    else
        SSP = nothing
    end

    # --------------------------------------------------------------------------    
    # Model Construction
    # --------------------------------------------------------------------------    

    # component first and lasts
    model_first = 1750
    brick_first = 1850
    damages_first = 2020
    model_last = 2300

    # Start with an instance of the FAIR model.
    m = MimiFAIRv1_6_2.get_model(start_year=model_first, end_year=model_last, ar6_scenario=ar6_scenario)

    # Set Dimensions
    set_dimension!(m, :time, model_first:model_last) # used in all components - already set in FAIR but reset for clarity
    set_dimension!(m, :country, countries) # used in most components

    set_dimension!(m, :fund_regions, fund_regions) # Agriculture components
    set_dimension!(m, :segments, segment_fingerprints.segments) # BRICK components
    set_dimension!(m, :ag_mapping_input_regions, countries) # Agriculture Aggregator components
    set_dimension!(m, :ag_mapping_output_regions, fund_regions) # Agriculture Aggregator components
    set_dimension!(m, :energy_countries, countries) # Countries used in energy damage function

    set_dimension!(m, :domestic_countries, domestic_countries) # Country ISO3 codes to be accumulated for domestic

    # Add Socioeconomics component BEFORE the FAIR model to allow for emissions feedbacks after damages_first year
    if socioeconomics_source == :RFF
        add_comp!(m, MimiRFFSPs.SPs, :Socioeconomic, first=damages_first, before=:ch4_cycle)
    elseif socioeconomics_source == :SSP
        add_comp!(m, MimiSSPs.SSPs, :Socioeconomic, first=damages_first, before=:ch4_cycle)
    end

    # Add PerCapitaGDP component
    add_comp!(m, PerCapitaGDP, :PerCapitaGDP, first=damages_first, after=:Socioeconomic)

    # Add VSL component
    add_comp!(m, VSL, :VSL, first=damages_first, after=:PerCapitaGDP)

    # We add an identity component that simply passes values through here
    # This makes it easier to later insert the marginal emission modification component
    # between two components that don't use backup data
    add_comp!(m, IdentityComponent_co2, :co2_emissions_identity, before=:co2_cycle)
    add_comp!(m, IdentityComponent_ch4, :ch4_emissions_identity, before=:ch4_cycle)
    add_comp!(m, IdentityComponent_n2o, :n2o_emissions_identity, before=:n2o_cycle)

    # Add Temperature Normalization Components
    add_comp!(m, GlobalTempNorm, :TempNorm_1880, after=:temperature) # Howard and Sterner
    add_comp!(m, GlobalTempNorm, :TempNorm_1900, after=:TempNorm_1880) # DICE
    add_comp!(m, GlobalTempNorm, :TempNorm_1850to1900, after=:TempNorm_1900) # Useful Reference to IPCC
    add_comp!(m, GlobalTempNorm, :TempNorm_1995to2005, after=:TempNorm_1850to1900) # Agriculture

    # Add Ocean Heat Accumulator to Link FAIR and BRICK
    add_comp!(m, OceanHeatAccumulator, after=:TempNorm_1995to2005)

    # Add BRICK components
    add_comp!(m, MimiBRICK.antarctic_ocean, first=brick_first, after=:OceanHeatAccumulator)
    add_comp!(m, MimiBRICK.antarctic_icesheet, first=brick_first, after=:antarctic_ocean)
    add_comp!(m, MimiBRICK.glaciers_small_icecaps, first=brick_first, after=:antarctic_icesheet)
    add_comp!(m, MimiBRICK.greenland_icesheet, first=brick_first, after=:glaciers_small_icecaps)
    add_comp!(m, MimiBRICK.thermal_expansion, first=brick_first, after=:greenland_icesheet)
    add_comp!(m, MimiBRICK.landwater_storage, first=brick_first, after=:thermal_expansion)
    add_comp!(m, MimiBRICK.global_sea_level, first=brick_first, after=:landwater_storage)

    # Add SLR Normalization components
    add_comp!(m, GlobalSLRNorm, :GlobalSLRNorm_1900, first=brick_first, after=:global_sea_level)

    # Add OceanPH components
    add_comp!(m, Mimi_NAS_pH.ocean_pH, :OceanPH, after=:GlobalSLRNorm_1900)

    # Add CromarMortality component
    add_comp!(m, cromar_mortality_damages, :CromarMortality, first=damages_first, after=:OceanPH)

    # Add Agriculture components
    add_comp!(m, Agriculture_RegionAggregatorSum, :Agriculture_aggregator_population, first=damages_first, after=:CromarMortality)
    add_comp!(m, Agriculture_RegionAggregatorSum, :Agriculture_aggregator_gdp, first=damages_first, after=:Agriculture_aggregator_population)
    add_comp!(m, MimiMooreEtAlAgricultureImpacts.Agriculture, :Agriculture, first=damages_first, after=:Agriculture_aggregator_gdp)
    add_comp!(m, AgricultureDamagesDisaggregator, :AgricultureDamagesDisaggregator, first=damages_first, after=:Agriculture)

    # add aggregators for 1990 population and GDP if we are using the RFFSPs
    socioeconomics_source == :RFF ? add_comp!(m, Agriculture_RegionAggregatorSum_NoTime, :Agriculture_aggregator_pop90, first=damages_first, after=:Agriculture_aggregator_gdp) : nothing
    socioeconomics_source == :RFF ? add_comp!(m, Agriculture_RegionAggregatorSum_NoTime, :Agriculture_aggregator_gdp90, first=damages_first, after=:Agriculture_aggregator_pop90) : nothing

    # Add a Regional Per Capita GDP component that takes inputs from the Agiculture aggregator components
    add_comp!(m, RegionalPerCapitaGDP, :RegionalPerCapitaGDP, first=damages_first, after=:AgricultureDamagesDisaggregator)

    # Add Energy components
    add_comp!(m, energy_damages, :energy_damages, first=damages_first, after=:RegionalPerCapitaGDP)

    # Add DICE2016R2 damage component
    add_comp!(m, dice2016R2_damage, :dice2016R2_damage, first=damages_first, after=:energy_damages)

    # Add Howard and Sterner damage components
    add_comp!(m, hs_damage, :hs_damage, first=damages_first, after=:dice2016R2_damage)

    # Add DamageAggregator component and regional damages aggregator helper function
    add_comp!(m, Damages_RegionAggregatorSum, first=damages_first)
    add_comp!(m, DamageAggregator, first=damages_first)

    # Add net consumption components (global and regional)
    add_comp!(m, GlobalNetConsumption, :global_netconsumption, first=damages_first, after=:DamageAggregator)
    add_comp!(m, RegionalNetConsumption, :regional_netconsumption, first=damages_first, after=:global_netconsumption)
    add_comp!(m, CountryNetConsumption, :country_netconsumption, first=damages_first, after=:regional_netconsumption)

    # --------------------------------------------------------------------------
    # Shared Model Parameters
    # --------------------------------------------------------------------------

    add_shared_param!(m, :model_country_names, countries, dims=[:country])
    add_shared_param!(m, :model_fund_regions, fund_regions, dims=[:fund_regions])

    # Agriculture
    add_shared_param!(m, :model_ag_mapping_input_regions, countries, dims=[:ag_mapping_input_regions])
    add_shared_param!(m, :model_ag_mapping_output_regions, fund_regions, dims=[:ag_mapping_output_regions])
    add_shared_param!(m, :model_ag_mapping, ag_mapping, dims=[:ag_mapping_input_regions])

    # BRICK
    add_shared_param!(m, :model_brick_seawater_freeze, -1.8)

    # Mortality
    if socioeconomics_source == :SSP
        mortality_data = load(joinpath(@__DIR__, "..", "data", "Mortality_cdr_spp_country_extensions_annual.csv")) |>
                         DataFrame |>
                         @filter(_.year in damages_first:model_last && _.scenario == mortality_SSP_map[SSP]) |>
                         DataFrame |>
                         @select(:year, :ISO, :cdf) |>
                         DataFrame |>
                         @orderby(:ISO) |>
                         DataFrame |>
                         i -> unstack(i, :year, :ISO, :cdf) |>
                              DataFrame |>
                              i -> select!(i, Not(:year))

        # make sure the columns match the mortality countries
        names(mortality_data) == countries ? nothing : "Countries in mortality data must match model countries."
        add_shared_param!(m, :model_ssp_baseline_mortality_rate, vcat(fill(NaN, (length(model_first:damages_first-1), size(mortality_data)[2])), mortality_data |> Matrix), dims=[:time, :country]) # Pad with NaN b/c starting component in later year.
    end

    # --------------------------------------------------------------------------
    # Component-Specific Parameters and Connections
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------    
    # BRICK
    # --------------------------------------------------------------------------

    # ----- Ocean Heat Accumulator ----- #

    connect_param!(m, :OceanHeatAccumulator, :del_ohc, :temperature, :del_ohc)

    # ----- Antarctic Ocean ----- #

    update_param!(m, :antarctic_ocean, :anto_α, 0.28)
    update_param!(m, :antarctic_ocean, :anto_β, 0.95)

    # ----- Antarctic Ice Sheet ----- #

    update_param!(m, :antarctic_icesheet, :ais_ρ_ice, 917.0)
    update_param!(m, :antarctic_icesheet, :ais_ρ_seawater, 1030.0)
    update_param!(m, :antarctic_icesheet, :ais_ρ_rock, 4000.0)
    update_param!(m, :antarctic_icesheet, :ais_sea_level₀, 0.0)
    update_param!(m, :antarctic_icesheet, :ais_ocean_temperature₀, 0.72)
    update_param!(m, :antarctic_icesheet, :ais_radius₀, 1.864e6)
    update_param!(m, :antarctic_icesheet, :ais_bedheight₀, 781.0)
    update_param!(m, :antarctic_icesheet, :ais_slope, 0.0006)
    update_param!(m, :antarctic_icesheet, :ais_μ, 11.0)
    update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, 1400.0)
    update_param!(m, :antarctic_icesheet, :ais_c, 100.0)
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀, 0.37)
    update_param!(m, :antarctic_icesheet, :ais_κ, 0.062)
    update_param!(m, :antarctic_icesheet, :ais_ν, 0.0086)
    update_param!(m, :antarctic_icesheet, :ais_iceflow₀, 1.2)
    update_param!(m, :antarctic_icesheet, :ais_γ, 2.9)
    update_param!(m, :antarctic_icesheet, :ais_α, 0.23)
    update_param!(m, :antarctic_icesheet, :ais_temperature_coefficient, 0.8365)
    update_param!(m, :antarctic_icesheet, :ais_temperature_intercept, 15.42)
    update_param!(m, :antarctic_icesheet, :ais_local_fingerprint, -1.18)
    update_param!(m, :antarctic_icesheet, :ocean_surface_area, 3.619e14)
    update_param!(m, :antarctic_icesheet, :temperature_threshold, -15.0)
    update_param!(m, :antarctic_icesheet, :λ, 0.0093)
    update_param!(m, :antarctic_icesheet, :include_ais_DSL, true)

    # ----- Glaciers & Small Ice Caps ----- #

    update_param!(m, :glaciers_small_icecaps, :gsic_β₀, 0.0013)
    update_param!(m, :glaciers_small_icecaps, :gsic_v₀, 0.376)
    update_param!(m, :glaciers_small_icecaps, :gsic_s₀, -0.0138)
    update_param!(m, :glaciers_small_icecaps, :gsic_n, 0.847)
    update_param!(m, :glaciers_small_icecaps, :gsic_teq, -0.15)

    # ----- Greenland Ice Sheet ----- #

    update_param!(m, :greenland_icesheet, :greenland_a, -1.37)
    update_param!(m, :greenland_icesheet, :greenland_b, 8.06)
    update_param!(m, :greenland_icesheet, :greenland_α, 0.0008)
    update_param!(m, :greenland_icesheet, :greenland_β, 0.00009)
    update_param!(m, :greenland_icesheet, :greenland_v₀, 7.52)

    # ----- Thermal Expansion ----- #

    update_param!(m, :thermal_expansion, :te_A, 3.619e14)
    update_param!(m, :thermal_expansion, :te_C, 3991.86795711963)
    update_param!(m, :thermal_expansion, :te_ρ, 1027.0)
    update_param!(m, :thermal_expansion, :te_α, 0.16)
    update_param!(m, :thermal_expansion, :te_s₀, 0.0)

    update_param!(m, :thermal_expansion, :ocean_heat_mixed, zeros(length(model_first:model_last)))
    connect_param!(m, :thermal_expansion, :ocean_heat_interior, :OceanHeatAccumulator, :del_ohc_accum)

    # ----- Landwater Storage ----- #

    update_param!(m, :landwater_storage, :lws₀, 0.0)
    update_param!(m, :landwater_storage, :first_projection_year, 2018)
    update_param!(m, :landwater_storage, :lws_random_sample, fill(0.0003, model_last - model_first + 1))

    # ----- Set Parameters With Common Values Across Components ----- #

    connect_param!(m, :antarctic_icesheet, :seawater_freeze, :model_brick_seawater_freeze)
    connect_param!(m, :antarctic_ocean, :seawater_freeze, :model_brick_seawater_freeze)

    update_param!(m, :GlobalSLRNorm_1900, :norm_range_start, 1900)
    update_param!(m, :GlobalSLRNorm_1900, :norm_range_end, 1900)

    # --------------------------------------------------------------------------    
    # Create Component Connections

    connect_param!(m, :global_sea_level => :slr_glaciers_small_ice_caps, :glaciers_small_icecaps => :gsic_sea_level)
    connect_param!(m, :global_sea_level => :slr_greeland_icesheet, :greenland_icesheet => :greenland_sea_level)
    connect_param!(m, :global_sea_level => :slr_antartic_icesheet, :antarctic_icesheet => :ais_sea_level)
    connect_param!(m, :global_sea_level => :slr_thermal_expansion, :thermal_expansion => :te_sea_level)
    connect_param!(m, :global_sea_level => :slr_landwater_storage, :landwater_storage => :lws_sea_level)

    connect_param!(m, :antarctic_icesheet => :antarctic_ocean_temperature, :antarctic_ocean => :anto_temperature)
    connect_param!(m, :antarctic_icesheet => :global_sea_level, :global_sea_level => :sea_level_rise)

    connect_param!(m, :antarctic_icesheet => :global_surface_temperature, :temperature => :T)
    connect_param!(m, :antarctic_ocean => :global_surface_temperature, :temperature => :T)
    connect_param!(m, :glaciers_small_icecaps => :global_surface_temperature, :temperature => :T)
    connect_param!(m, :greenland_icesheet => :global_surface_temperature, :temperature => :T)

    connect_param!(m, :GlobalSLRNorm_1900 => :global_slr, :global_sea_level => :sea_level_rise)

    # --------------------------------------------------------------------------    
    # OceanPH
    # --------------------------------------------------------------------------

    update_param!(m, :OceanPH, :β1, -0.3671)
    update_param!(m, :OceanPH, :β2, 10.2328)
    update_param!(m, :OceanPH, :pH_0, 8.123)

    connect_param!(m, :OceanPH => :atm_co2_conc, :co2_cycle => :co2)

    # --------------------------------------------------------------------------    
    # Socioeconomic
    # --------------------------------------------------------------------------

    if socioeconomics_source == :SSP
        update_param!(m, :Socioeconomic, :SSP_source, "Benveniste") # only available source to 2300 at this time in MimiSSPs
        update_param!(m, :Socioeconomic, :SSP, SSP) # select the SSP from RCMIP name ie. SSP2
        update_param!(m, :Socioeconomic, :emissions_source, "Leach") # only available source to 2300 at this time in MimiSSPs
        update_param!(m, :Socioeconomic, :emissions_scenario, SSP_scenario) # full name ie. SSSP245

    elseif socioeconomics_source == :RFF
        isnothing(RFFSPsample) ? nothing : update_param!(m, :Socioeconomic, :id, RFFSPsample)
    end
    connect_param!(m, :Socioeconomic, :country_names, :model_country_names)

    # Feedback of Socioeconomic Emissions back to FAIR

    # Load IPCC AR6 emissions scenario used for FAIRv1.6.2 ensemble runs (options = "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585").
    ar6_emissions_raw = DataFrame(load(joinpath(@__DIR__, "..", "data", "FAIR_ar6", "AR6_emissions_" * ar6_scenario * "_1750_2300.csv")))

    # Subset AR6 emissions to proper years.
    emission_indices = indexin(collect(model_first:model_last), ar6_emissions_raw.Year)
    ar6_emissions = ar6_emissions_raw[emission_indices, :]

    # Here we couple the identity component co2_emissions to the SSP output, and then the
    # FAIR emissions component to that identity component co2_emissions
    connect_param!(m, :co2_emissions_identity => :input_co2, :Socioeconomic => :co2_emissions, ar6_emissions.FossilCO2 .+ ar6_emissions.OtherCO2)
    connect_param!(m, :co2_cycle => :E_co2, :co2_emissions_identity => :output_co2)

    # do the same for n2o_emissions
    connect_param!(m, :n2o_emissions_identity => :input_n2o, :Socioeconomic => :n2o_emissions, ar6_emissions.N2O)
    connect_param!(m, :n2o_cycle => :fossil_emiss_N₂O, :n2o_emissions_identity => :output_n2o)

    # do the same for ch4_emissions
    connect_param!(m, :ch4_emissions_identity => :input_ch4, :Socioeconomic => :ch4_emissions, ar6_emissions.CH4)
    connect_param!(m, :ch4_cycle => :fossil_emiss_CH₄, :ch4_emissions_identity => :output_ch4)

    # Land Use CO2 Emissions - FAIRv1.6.2 component :landuse_forcing and parameter :landuse_emiss
    #
    # For the SSPs, the MimiSSPs component does not break carbon dioxide out by industrial 
    # and other, so we will simply let FAIR1.6.2 run with its original settings for land use CO2, which is 
    # consistent with Leach (FAIRv2.0) but just not broken out in that dataset, so this is consistent. 
    # For the RFF SPs we can either let them run with the middle-of the road and best matched (pop and emissions)
    # ssp245 AR6 scenario, or explicitly connect new data. Currently do the former.

    # --------------------------------------------------------------------------    
    # PerCapitaGDP
    # -------------------------------------------------------------------------- 

    connect_param!(m, :PerCapitaGDP => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :PerCapitaGDP => :population, :Socioeconomic => :population)

    # --------------------------------------------------------------------------    
    # VSL
    # --------------------------------------------------------------------------

    if vsl == :fund
        update_param!(m, :VSL, :α, 4.99252262888626e6 * pricelevel_1995_to_2005)   # convert from FUND USD $1995 to USD $2005
        update_param!(m, :VSL, :y₀, 24_962.6131444313 * pricelevel_1995_to_2005)   # convert from FUND USD $1995 to USD $2005
    elseif vsl == :epa
        update_param!(m, :VSL, :α, 7.73514707e6)                                   # 2020 EPA VSL in 2005$. See DataExplainer.ipynb for information
        update_param!(m, :VSL, :y₀, 48_726.60)                                      # 2020 U.S. income per capita in 2005$; See DataExplainer.ipynb for information  
    elseif vsl == :uba
        update_param!(m, :VSL, :α, 5_920_000. / pricelevel_2005_to_2020)           # 2020 UBA VSL in 2005$
        update_param!(m, :VSL, :y₀, 44_646.78)                                      # 2020 German income per capita in 2005$
    else
        error("Invalid vsl argument of $vsl.")
    end

    update_param!(m, :VSL, :ϵ, 1.0)
    connect_param!(m, :VSL => :pc_gdp, :PerCapitaGDP => :pc_gdp)

    # --------------------------------------------------------------------------    
    # Temperature Normlization Components
    # --------------------------------------------------------------------------

    # TempNorm_1880 - Normalize temperature to deviation from 1880 for Howard and Sterner damage function
    update_param!(m, :TempNorm_1880, :norm_range_start, 1880)
    update_param!(m, :TempNorm_1880, :norm_range_end, 1880)
    connect_param!(m, :TempNorm_1880 => :global_temperature, :temperature => :T)

    # TempNorm_1900 - Normalize temperature to deviation from 1900 for DICE2016 damage function
    update_param!(m, :TempNorm_1900, :norm_range_start, 1900)
    update_param!(m, :TempNorm_1900, :norm_range_end, 1900)
    connect_param!(m, :TempNorm_1900 => :global_temperature, :temperature => :T)

    # TempNorm_1850to1900 - Normalize temperature to deviation from 1850 to 1900 for IPCC Comparison Graphics
    update_param!(m, :TempNorm_1850to1900, :norm_range_start, 1850)
    update_param!(m, :TempNorm_1850to1900, :norm_range_end, 1900)
    connect_param!(m, :TempNorm_1850to1900 => :global_temperature, :temperature => :T)

    # TempNorm_1995to2005 - Normalize temperature to deviation from 1995 to 2005 for Agriculture Component
    update_param!(m, :TempNorm_1995to2005, :norm_range_start, 1995)
    update_param!(m, :TempNorm_1995to2005, :norm_range_end, 2005)
    connect_param!(m, :TempNorm_1995to2005 => :global_temperature, :temperature => :T)

    # --------------------------------------------------------------------------
    # Cromar et al. Temperature-Mortality Damages
    # --------------------------------------------------------------------------

    # Assign Cromar et al. regional temperature mortality coefficients to appropriate countries.

    # Load raw data.
    cromar_coeffs = load(joinpath(@__DIR__, "..", "data", "CromarMortality_damages_coefficients.csv")) |> DataFrame
    cromar_mapping_raw = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_cromar_mortality_regions.csv")) |> DataFrame

    # Initialize an array to store country-level coefficients
    country_β_mortality = zeros(length(cromar_mapping_raw.ISO3))

    # Loop through the regions and assign regional coefficients to proper sets of countries.
    for r = 1:length(cromar_regions)
        # Find country indices for region "r"
        r_index = findall(x -> x == cromar_regions[r], cromar_mapping_raw.cromar_region)
        # Find index for region "r" coefficient.
        β_index = findfirst(x -> x == cromar_regions[r], cromar_coeffs[!, "Cromar Region Name"])
        # Assign all countries in that region proper coefficient.
        country_β_mortality[r_index] .= cromar_coeffs[β_index, "Pooled Beta"]
    end

    # Get indices to reorder Cromar countries mapped to countries dimension (could be correct oder already, this is a safety check)
    cromar_indices = indexin(countries, cromar_mapping_raw.ISO3)
    country_β_mortality = country_β_mortality[cromar_indices]

    update_param!(m, :CromarMortality, :β_mortality, country_β_mortality)

    if socioeconomics_source == :SSP
        connect_param!(m, :CromarMortality, :baseline_mortality_rate, :model_ssp_baseline_mortality_rate) # shared model parameter
    elseif socioeconomics_source == :RFF
        connect_param!(m, :CromarMortality => :baseline_mortality_rate, :Socioeconomic => :deathrate)
    end

    connect_param!(m, :CromarMortality => :population, :Socioeconomic => :population)
    connect_param!(m, :CromarMortality => :temperature, :temperature => :T)
    connect_param!(m, :CromarMortality => :vsl, :VSL => :vsl)

    # --------------------------------------------------------------------------
    # Agriculture Aggregators
    # --------------------------------------------------------------------------

    connect_param!(m, :Agriculture_aggregator_population, :input_region_names, :model_ag_mapping_input_regions)
    connect_param!(m, :Agriculture_aggregator_population, :output_region_names, :model_ag_mapping_output_regions)
    connect_param!(m, :Agriculture_aggregator_population, :input_output_mapping, :model_ag_mapping)
    connect_param!(m, :Agriculture_aggregator_population => :input, :Socioeconomic => :population)

    connect_param!(m, :Agriculture_aggregator_gdp, :input_region_names, :model_ag_mapping_input_regions)
    connect_param!(m, :Agriculture_aggregator_gdp, :output_region_names, :model_ag_mapping_output_regions)
    connect_param!(m, :Agriculture_aggregator_gdp, :input_output_mapping, :model_ag_mapping)
    connect_param!(m, :Agriculture_aggregator_gdp => :input, :Socioeconomic => :gdp)

    if socioeconomics_source == :RFF
        connect_param!(m, :Agriculture_aggregator_gdp90, :input_region_names, :model_ag_mapping_input_regions)
        connect_param!(m, :Agriculture_aggregator_gdp90, :output_region_names, :model_ag_mapping_output_regions)
        connect_param!(m, :Agriculture_aggregator_gdp90, :input_output_mapping, :model_ag_mapping)
        connect_param!(m, :Agriculture_aggregator_gdp90 => :input, :Socioeconomic => :gdp1990)

        connect_param!(m, :Agriculture_aggregator_pop90, :input_region_names, :model_ag_mapping_input_regions)
        connect_param!(m, :Agriculture_aggregator_pop90, :output_region_names, :model_ag_mapping_output_regions)
        connect_param!(m, :Agriculture_aggregator_pop90, :input_output_mapping, :model_ag_mapping)
        connect_param!(m, :Agriculture_aggregator_pop90 => :input, :Socioeconomic => :population1990)
    end

    # --------------------------------------------------------------------------
    # Agriculture
    # --------------------------------------------------------------------------

    fund_regions != MimiMooreEtAlAgricultureImpacts.fund_regions && error("FUND regions for RFF Model do not match FUND regions for Agriculture.")

    # Handle in pop and gdp 1990 baseline values
    if socioeconomics_source == :SSP

        data1990 = load(joinpath(@__DIR__, "..", "data", "Benveniste_SSPs", "Agriculture_1990vals.csv")) |>
                   DataFrame |>
                   @filter(_.SSP == SSP) |>
                   DataFrame
        idxs = indexin(data1990.fund_region, fund_regions) # get the ordering of 1990 regions matched to fund regions in model
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("FUND regions for RFF Model do not match FUND regions for Agriculture 1990 values.") : nothing
        data1990 = data1990[idxs, :] # reorder based on idxs

        update_param!(m, :Agriculture, :pop90, data1990.pop)
        update_param!(m, :Agriculture, :gdp90, data1990.gdp)

    elseif socioeconomics_source == :RFF
        connect_param!(m, :Agriculture => :pop90, :Agriculture_aggregator_pop90 => :output)
        connect_param!(m, :Agriculture => :gdp90, :Agriculture_aggregator_gdp90 => :output)
    end

    # Access which of the 5 possible DFs to use for the damage function
    gtap_idx = findfirst(isequal(Agriculture_gtap), MimiMooreEtAlAgricultureImpacts.gtaps)
    gtap_df = MimiMooreEtAlAgricultureImpacts.gtap_df_all[:, :, gtap_idx]

    update_param!(m, :Agriculture, :gtap_df, gtap_df)
    update_param!(m, :Agriculture, :gtap_name, Agriculture_gtap)
    update_param!(m, :Agriculture, :floor_on_damages, Agriculture_floor_on_damages)
    update_param!(m, :Agriculture, :ceiling_on_benefits, Agriculture_ceiling_on_benefits)
    update_param!(m, :Agriculture, :agrish0, Array{Float64,1}(readdlm(joinpath(MimiMooreEtAlAgricultureImpacts.fund_datadir, "agrish0.csv"), ',', skipstart=1)[:, 2]))

    connect_param!(m, :Agriculture => :population, :Agriculture_aggregator_population => :output)
    connect_param!(m, :Agriculture => :income, :Agriculture_aggregator_gdp => :output)
    connect_param!(m, :Agriculture => :temp, :TempNorm_1995to2005 => :global_temperature_norm)

    # --------------------------------------------------------------------------
    # Regional Per Capita GDP
    # --------------------------------------------------------------------------

    connect_param!(m, :RegionalPerCapitaGDP => :population, :Agriculture_aggregator_population => :output)
    connect_param!(m, :RegionalPerCapitaGDP => :gdp, :Agriculture_aggregator_gdp => :output)

    # --------------------------------------------------------------------------
    # Agriculture Damages Disaggregator
    # --------------------------------------------------------------------------

    connect_param!(m, :AgricultureDamagesDisaggregator, :mapping, :model_ag_mapping)
    connect_param!(m, :AgricultureDamagesDisaggregator, :fund_region_names, :model_ag_mapping_output_regions)

    connect_param!(m, :AgricultureDamagesDisaggregator => :gdp_fund_region, :Agriculture_aggregator_gdp => :output)
    connect_param!(m, :AgricultureDamagesDisaggregator => :gdp_country, :Socioeconomic => :gdp)

    connect_param!(m, :AgricultureDamagesDisaggregator => :damages_ag_fund_region, :Agriculture => :agcost)

    # --------------------------------------------------------------------------
    # Energy
    # --------------------------------------------------------------------------

    # Assign GCAM regional energy damage coefficients to appropriate countries.

    # Load raw data.
    energy_coeffs = load(joinpath(@__DIR__, "..", "data", "energy_damages_gcam_region_coefficients.csv")) |> DataFrame
    gcam_mapping_raw = load(joinpath(@__DIR__, "..", "data", "Mapping_countries_to_gcam_energy_regions.csv")) |> DataFrame

    # Initialize an array to store country-level coefficients
    country_β_energy = zeros(length(gcam_mapping_raw.ISO3))

    # Loop through the regions and assign regional coefficients to proper subset of countries.
    for r = 1:length(gcam_regions)
        # Find country indices for region "r"
        r_index = findall(x -> x == gcam_regions[r], gcam_mapping_raw.gcamregion)
        # Find index for region "r" coefficient.
        β_index = findfirst(x -> x == gcam_regions[r], energy_coeffs.gcam_region)
        # Assign all countries in that region proper coefficient.
        country_β_energy[r_index] .= energy_coeffs[β_index, "coefficient"]
    end

    set_param!(m, :energy_damages, :β_energy, country_β_energy)
    connect_param!(m, :energy_damages => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :energy_damages => :temperature, :temperature => :T)

    # --------------------------------------------------------------------------
    # DICE2016R2 Damages
    # --------------------------------------------------------------------------

    connect_param!(m, :dice2016R2_damage => :temperature, :TempNorm_1900 => :global_temperature_norm)
    connect_param!(m, :dice2016R2_damage => :gdp, :Socioeconomic => :gdp)

    # --------------------------------------------------------------------------
    # Howard and Sterner Damages
    # --------------------------------------------------------------------------

    connect_param!(m, :hs_damage => :temperature, :TempNorm_1880 => :global_temperature_norm)
    connect_param!(m, :hs_damage => :gdp, :Socioeconomic => :gdp)

    # --------------------------------------------------------------------------
    # Damage Aggregation
    # --------------------------------------------------------------------------

    # small regional damage aggregator helper component
    connect_param!(m, :Damages_RegionAggregatorSum, :input_region_names, :model_ag_mapping_input_regions)
    connect_param!(m, :Damages_RegionAggregatorSum, :output_region_names, :model_ag_mapping_output_regions)
    connect_param!(m, :Damages_RegionAggregatorSum, :input_output_mapping, :model_ag_mapping)
    connect_param!(m, :Damages_RegionAggregatorSum => :damage_cromar_mortality, :CromarMortality => :mortality_costs)
    connect_param!(m, :Damages_RegionAggregatorSum => :damage_energy, :energy_damages => :energy_costs_dollar)

    # main damage aggregator
    connect_param!(m, :DamageAggregator => :damage_ag, :Agriculture => :agcost)
    connect_param!(m, :DamageAggregator => :damage_ag_countries, :AgricultureDamagesDisaggregator => :damages_ag_country)
    connect_param!(m, :DamageAggregator => :damage_cromar_mortality, :CromarMortality => :mortality_costs)
    connect_param!(m, :DamageAggregator => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :DamageAggregator => :damage_energy, :energy_damages => :energy_costs_dollar)
    connect_param!(m, :DamageAggregator => :damage_dice2016R2, :dice2016R2_damage => :damages)
    connect_param!(m, :DamageAggregator => :damage_hs, :hs_damage => :damages)

    connect_param!(m, :DamageAggregator => :damage_cromar_mortality_regions, :Damages_RegionAggregatorSum => :damage_cromar_mortality_regions)
    connect_param!(m, :DamageAggregator => :damage_energy_regions, :Damages_RegionAggregatorSum => :damage_energy_regions)

    domestic_idxs_country_dim = Int.(indexin(dim_keys(m, :domestic_countries), dim_keys(m, :country)))
    update_param!(m, :DamageAggregator, :domestic_idxs_country_dim, domestic_idxs_country_dim)

    domestic_idxs_energy_countries_dim = Int.(indexin(dim_keys(m, :domestic_countries), dim_keys(m, :energy_countries)))
    update_param!(m, :DamageAggregator, :domestic_idxs_energy_countries_dim, domestic_idxs_energy_countries_dim)

    # --------------------------------------------------------------------------
    # Net Consumption
    # --------------------------------------------------------------------------

    # global
    connect_param!(m, :global_netconsumption => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :global_netconsumption => :population, :Socioeconomic => :population)
    connect_param!(m, :global_netconsumption => :total_damage, :DamageAggregator => :total_damage)

    # regional
    connect_param!(m, :regional_netconsumption => :population, :Agriculture_aggregator_population => :output)
    connect_param!(m, :regional_netconsumption => :gdp, :Agriculture_aggregator_gdp => :output)
    connect_param!(m, :regional_netconsumption => :total_damage, :DamageAggregator => :total_damage_regions)

    # country
    connect_param!(m, :country_netconsumption => :gdp, :Socioeconomic => :gdp)
    connect_param!(m, :country_netconsumption => :population, :Socioeconomic => :population)
    connect_param!(m, :country_netconsumption => :total_damage, :DamageAggregator => :total_damage_countries)

    return m
end
