using Mimi

# ------------------------------------------------------------------------------
# Calculate the Per Capita GDP
# ------------------------------------------------------------------------------

@defcomp PerCapitaGDP begin

    country = Index()

    gdp         = Parameter(index=[time, country], unit="billion US\$2005/yr")
    population  = Parameter(index=[time, country], unit="million")

    pc_gdp      = Variable(index=[time, country], unit = "US\$2005/yr/person") # Country-level per capita GDP ($/person).
    global_pc_gdp       = Variable(index=[time], unit = "US\$2005/yr/person")

    function run_timestep(p, v, d, t)

        # calculate global per capita gdp
        v.global_pc_gdp[t] = sum(p.gdp[t,:]) / sum(p.population[t,:]) * 1e3 

        # calculate country level per capita gdp
        for c in d.country
            v.pc_gdp[t, c] = (p.gdp[t, c]) ./ (p.population[t, c]) * 1e3 
        end
    end
end
