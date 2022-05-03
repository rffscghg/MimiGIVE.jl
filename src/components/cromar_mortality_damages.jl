using Mimi

# ------------------------------------------------------------
# Temperature Mortality Damages (based on Cromar et al. 2021)
# ------------------------------------------------------------

@defcomp cromar_mortality_damages begin

    country                 = Index() # Index for countries in the regions used for the Cromar et al. temperature-mortality damage functions.

   	β_mortality             = Parameter(index=[country]) # Coefficient relating global temperature to change in mortality rates.
    baseline_mortality_rate = Parameter(index=[time, country]) # Crude death rate in a given country (deaths per 1,000 population).
 	temperature             = Parameter(index=[time], unit="degC") # Global average surface temperature anomaly relative to pre-industrial (°C).

    population              = Parameter(index=[time, country], unit="million") # Population in a given country (millions of persons).
    vsl                     = Parameter(index=[time, country], unit="US\$2005/yr") # Value of a statistical life ($).

    mortality_change        = Variable(index=[time, country])  # Change in a country's baseline mortality rate due to combined effects of cold and heat (with positive values indicating increasing mortality rates).
   	mortality_costs         = Variable(index=[time, country], unit="US\$2005/yr")  # Costs of temperature mortality based on the VSL ($).
    excess_death_rate       = Variable(index=[time, country], unit = "deaths/1000 persons/yr")  # Change in a country's baseline death rate due to combined effects of cold and heat (additional deaths per 1,000 population).
    excess_deaths           = Variable(index=[time, country], unit="persons")  # Additional deaths that occur in a country due to the combined effects of cold and heat (individual persons).


      function run_timestep(p, v, d, t)

        for c in d.country

            # Calculate change in a country's baseline mortality rate due to combined effects of heat and cold.
            v.mortality_change[t,c] = p.β_mortality[c] * p.temperature[t]

            # Calculate additional deaths per 1,000 population due to cold and heat.
            v.excess_death_rate[t,c] = p.baseline_mortality_rate[t,c] * v.mortality_change[t,c]

            # Calculate additional deaths that occur due to cold and heat (assumes population units in millions of persons so converts to thousands to match deathrates units).
            v.excess_deaths[t,c] = (p.population[t,c] .* 1000) * v.excess_death_rate[t,c]

            # Multiply excess deaths by the VSL.
            v.mortality_costs[t,c] = p.vsl[t,c] * v.excess_deaths[t,c]
        end
    end
end
