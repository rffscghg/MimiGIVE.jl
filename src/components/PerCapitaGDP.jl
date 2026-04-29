using Mimi

# Calculate the per capita GDP

@defcomp PerCapitaGDP begin

    country = Index()

    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")
    population = Parameter(index=[time, country], unit="million")

    pc_gdp = Variable(index=[time, country], unit="US\$2005/yr/person") # Country-level per capita GDP ($/person).
    global_pc_gdp = Variable(index=[time], unit="US\$2005/yr/person")

    function run_timestep(p, v, d, t)

        # calculate global per capita gdp
        v.global_pc_gdp[t] = sum(p.gdp[t, :]) / sum(p.population[t, :]) * 1e3

        # calculate country level per capita gdp
        for c in d.country
            v.pc_gdp[t, c] = (p.gdp[t, c]) ./ (p.population[t, c]) * 1e3
        end
    end
end

@defcomp RegionalPerCapitaGDP begin

    fund_regions = Index()

    gdp = Parameter(index=[time, fund_regions], unit="billion US\$2005/yr")
    population = Parameter(index=[time, fund_regions], unit="million")

    pc_gdp = Variable(index=[time, fund_regions], unit="US\$2005/yr/person") # Region-level per capita GDP ($/person).

    function run_timestep(p, v, d, t)

        # calculate region level per capita gdp
        for r in d.fund_regions
            v.pc_gdp[t, r] = (p.gdp[t, r]) ./ (p.population[t, r]) * 1e3
        end
    end
end
