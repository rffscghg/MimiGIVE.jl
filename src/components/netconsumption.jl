using Mimi

@defcomp GlobalNetConsumption begin
    country = Index()

    gdp = Parameter(index=[time,country], unit="billion US\$2005/yr")
    population  = Parameter(index=[time, country], unit="million")
    total_damage = Parameter(index=[time], unit="US\$2005/yr")

    net_consumption = Variable(index=[time])
    net_cpc = Variable(index=[time])

    global_gdp = Variable(index=[time])
    global_population = Variable(index=[time])

    function run_timestep(p, v, d, t)

        # Sum the population and gdp of all countries for the current timestep
        v.global_population[t] = sum(p.population[t,:])
        v.global_gdp[t] = sum(p.gdp[t,:])

        # Convert damages to billions
        total_damage = p.total_damage[t] / 1e9

        # Compute net consumption as GDP - damages
        v.net_consumption[t] = v.global_gdp[t] - total_damage

        # We multiply by 1e3 because net_consumption is in billion, and population is in million
        v.net_cpc[t] = v.net_consumption[t] * 1e3 / v.global_population[t]

    end
end
