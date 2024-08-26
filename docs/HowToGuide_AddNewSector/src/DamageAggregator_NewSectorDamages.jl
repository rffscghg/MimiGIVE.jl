using Mimi

# Aggregate damages across damage functions

@defcomp DamageAggregator_NewSectorDamages begin

    fund_regions = Index()
    country = Index()
    energy_countries = Index()
    domestic_countries = Index()

    domestic_idxs_country_dim = Parameter{Int}(index=[domestic_countries])
    domestic_idxs_energy_countries_dim = Parameter{Int}(index=[domestic_countries])

    # internally compute for speed
    domestic_idxs_country_dim_int = Variable{Int}(index=[domestic_countries])
    domestic_idxs_energy_countries_dim_int = Variable{Int}(index=[domestic_countries])

    # inclusion of different damages

    # By default the individual sectoral damage calculations are ON, including 
    # SLR which runs after the main model, while global damage function calculations
    # are OFF.
    include_cromar_mortality = Parameter{Bool}(default=true)
    include_ag = Parameter{Bool}(default=true)
    include_slr = Parameter{Bool}(default=true)
    include_energy = Parameter{Bool}(default=true)
    include_new_sector = Parameter{Bool}(default=true)
    include_dice2016R2 = Parameter{Bool}(default=false)
    include_hs_damage = Parameter{Bool}(default=false)

    damage_cromar_mortality = Parameter(index=[time, country], unit="US\$2005/yr")
    damage_ag = Parameter(index=[time, fund_regions], unit="billion US\$2005/yr")
    damage_ag_countries = Parameter(index=[time, country], unit="billion US\$2005/yr") # ag damages disaggregated via method in AgricultureDamagesDisaggregator
    damage_energy = Parameter(index=[time, energy_countries], unit="billion US\$2005/yr")
    damage_new_sector = Parameter(index=[time, country], unit="billion US\$2005/yr")
    damage_dice2016R2 = Parameter(index=[time], unit="billion US\$2005/yr")
    damage_hs = Parameter(index=[time], unit="billion US\$2005/yr")

    # damages aggregated by fund regions
    damage_cromar_mortality_regions = Parameter(index=[time, fund_regions], unit="US\$2005/yr")
    damage_energy_regions = Parameter(index=[time, fund_regions], unit="billion US\$2005/yr")
    damage_new_sector_regions = Parameter(index=[time, fund_regions], unit="billion US\$2005/yr")

    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")

    total_damage = Variable(index=[time], unit="US\$2005/yr")
    total_damage_regions = Variable(index=[time, fund_regions], unit="US\$2005/yr")
    total_damage_countries = Variable(index=[time, country], unit="US\$2005/yr") # ag damages disaggregated via method in AgricultureDamagesDisaggregator
    total_damage_share = Variable(index=[time])
    total_damage_domestic = Variable(index=[time], unit="US\$2005/yr")

    # global annual aggregates - for interim model outputs and partial SCCs
    cromar_mortality_damage = Variable(index=[time], unit="US\$2005/yr")
    agriculture_damage = Variable(index=[time], unit="US\$2005/yr")
    energy_damage = Variable(index=[time], unit="US\$2005/yr")
    new_sector_damage = Variable(index=[time], unit="US\$2005/yr")

    # domestic annual aggregates - for interim model outputs and partial SCCs
    cromar_mortality_damage_domestic = Variable(index=[time], unit="US\$2005/yr")
    agriculture_damage_domestic = Variable(index=[time], unit="US\$2005/yr")
    energy_damage_domestic = Variable(index=[time], unit="US\$2005/yr")
    new_sector_damage_domestic = Variable(index=[time], unit="US\$2005/yr")

    function init(p, v, d)
        # convert to integers for indexing - do once here for speed
        v.domestic_idxs_country_dim_int[:] = Int.(p.domestic_idxs_country_dim)
        v.domestic_idxs_energy_countries_dim_int[:] = Int.(p.domestic_idxs_energy_countries_dim)
    end

    function run_timestep(p, v, d, t)

        # regional annual aggregates
        for r in d.fund_regions
            v.total_damage_regions[t, r] =
                (p.include_cromar_mortality ? p.damage_cromar_mortality_regions[t, r] : 0.) +
                (p.include_ag ? p.damage_ag[t, r] * 1e9 : 0.) +
                (p.include_energy ? p.damage_energy_regions[t, r] * 1e9 : 0.) +
                (p.include_new_sector ? p.damage_new_sector_regions[t, r] * 1e9 : 0.)

        end

        # country level aggregates where ag damages disaggregated via method in
        # AgricultureDamagesDisaggregator
        num_countries = length(d.country)
        v.total_damage_countries[t, :] =
            (p.include_cromar_mortality ? p.damage_cromar_mortality[t, :] : fill(0., num_countries)) +
            (p.include_ag ? p.damage_ag_countries[t, :] * 1e9 : fill(0., num_countries)) +
            (p.include_energy ? p.damage_energy[t, :] * 1e9 : fill(0., num_countries)) +
            (p.include_new_sector ? p.damage_new_sector[t, :] * 1e9 : fill(0., num_countries))

        # global annual aggregates - for interim model outputs and partial SCCs
        v.cromar_mortality_damage[t] = sum(p.damage_cromar_mortality[t, :])
        v.agriculture_damage[t] = sum(p.damage_ag[t, :]) * 1e9
        v.energy_damage[t] = sum(p.damage_energy[t, :]) * 1e9
        v.new_sector_damage[t] = sum(p.damage_new_sector[t, :]) * 1e9

        v.total_damage[t] =
            (p.include_cromar_mortality ? v.cromar_mortality_damage[t] : 0.) +
            (p.include_ag ? v.agriculture_damage[t] : 0.) +
            (p.include_energy ? v.energy_damage[t] : 0.) +
            (p.include_new_sector ? v.new_sector_damage[t] : 0.) +
            (p.include_dice2016R2 ? p.damage_dice2016R2[t] * 1e9 : 0.) +
            (p.include_hs_damage ? p.damage_hs[t] * 1e9 : 0.)

        gdp = sum(p.gdp[t, :]) * 1e9

        v.total_damage_share[t] = v.total_damage[t] / gdp

        # domestic annual aggregates - for interim model outputs and partial SCCs
        v.cromar_mortality_damage_domestic[t] = sum(p.damage_cromar_mortality[t, v.domestic_idxs_country_dim_int])
        v.agriculture_damage_domestic[t] = p.damage_ag[t, 1] * 1e9
        v.energy_damage_domestic[t] = sum(p.damage_energy[t, v.domestic_idxs_energy_countries_dim_int] * 1e9)
        v.new_sector_damage_domestic[t] = sum(p.damage_new_sector[t, v.domestic_idxs_country_dim_int] * 1e9)

        # Calculate domestic damages
        v.total_damage_domestic[t] =
            (p.include_cromar_mortality ? v.cromar_mortality_damage_domestic[t] : 0.) +
            (p.include_ag ? v.agriculture_damage_domestic[t] : 0.) +
            (p.include_energy ? v.energy_damage_domestic[t] : 0.) +
            (p.include_new_sector ? v.new_sector_damage_domestic[t] : 0.)
    end
end
