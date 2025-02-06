using Mimi

# Calculate global net consumption

@defcomp GlobalNetConsumption begin
    country = Index()

    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")
    population = Parameter(index=[time, country], unit="million")
    total_damage = Parameter(index=[time], unit="US\$2005/yr")

    net_consumption = Variable(index=[time])
    net_cpc = Variable(index=[time])

    global_gdp = Variable(index=[time])
    global_population = Variable(index=[time])

    function run_timestep(p, v, d, t)

        # Sum the population and gdp of all countries for the current timestep
        v.global_population[t] = sum(p.population[t, :])
        v.global_gdp[t] = sum(p.gdp[t, :])

        # Convert damages to billions
        total_damage = p.total_damage[t] / 1e9

        # Compute net consumption as GDP - damages
        v.net_consumption[t] = v.global_gdp[t] - total_damage

        # Multiply by 1e3 because net_consumption is in billion, and population is in million
        v.net_cpc[t] = v.net_consumption[t] * 1e3 / v.global_population[t]

    end
end

# Calculate regional net consumption

@defcomp RegionalNetConsumption begin
    fund_regions = Index()

    gdp = Parameter(index=[time, fund_regions], unit="billion US\$2005/yr")
    population = Parameter(index=[time, fund_regions], unit="million")
    total_damage = Parameter(index=[time, fund_regions], unit="US\$2005/yr")

    net_consumption = Variable(index=[time, fund_regions])
    net_cpc = Variable(index=[time, fund_regions])

    function run_timestep(p, v, d, t)

        for r in d.fund_regions
            # Convert damages to billions
            total_damage = p.total_damage[t, r] / 1e9

            # Compute net consumption as GDP - damages
            v.net_consumption[t, r] = p.gdp[t, r] - total_damage

            # Multiply by 1e3 because net_consumption is in billion, and population is in million
            v.net_cpc[t, r] = v.net_consumption[t, r] * 1e3 / p.population[t, r]
        end
    end
end

# Calculate country level net consumption

@defcomp CountryNetConsumption begin
    fund_regions = Index()

    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")
    population = Parameter(index=[time, country], unit="million")
    total_damage = Parameter(index=[time, country], unit="US\$2005/yr")

    net_consumption = Variable(index=[time, country])
    net_cpc = Variable(index=[time, country])

    function run_timestep(p, v, d, t)

        for c in d.country
            # Convert damages to billions
            total_damage = p.total_damage[t, c] / 1e9

            # Compute net consumption as GDP - damages
            v.net_consumption[t, c] = p.gdp[t, c] - total_damage

            # Multiply by 1e3 because net_consumption is in billion, and population is in million
            v.net_cpc[t, c] = v.net_consumption[t, c] * 1e3 / p.population[t, c]
        end
    end
end
