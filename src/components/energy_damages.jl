using Mimi

# --------------------------------------------------
# Energy-Use Damages (based on Clarke et al. 2018)
# --------------------------------------------------

@defcomp energy_damages begin

	energy_countries    = Index() 							        # Index for countries in the GCAM regions used for energy damage functions.

   	β_energy            = Parameter(index=[energy_countries])       # Coefficient relating global tempeature to change in energy expenditures as a share of GDP.
  	gdp 				= Parameter(index=[time, energy_countries], unit="billion US\$2005/yr") # Country-level GDP (billions US $2005 / yr").
 	temperature         = Parameter(index=[time], unit="degC") # Global average surface temperature anomaly relative to pre-industrial (°C).

   	energy_costs_dollar = Variable(index=[time, energy_countries], unit="billion US\$2005/yr")  # Change in energy expenditures in dollars (billions US $2005 / yr).
   	energy_costs_share  = Variable(index=[time, energy_countries])  # Change in energy expenditures as a share of GDP (Δ gdp share / °C).


    function run_timestep(p, v, d, t)

        for c in d.energy_countries

        	# Calculate additional energy expenditures as a share of GDP (coefficient gives percentages, so divide by 100 to get share).
        	v.energy_costs_share[t,c] = p.β_energy[c] * p.temperature[t] / 100.0

        	# Calculate additoinal energy expenditures in monetary terms.
        	v.energy_costs_dollar[t,c] = v.energy_costs_share[t,c] * p.gdp[t,c]
        end
    end
end
