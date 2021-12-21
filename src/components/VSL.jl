using Mimi

# ------------------------------------------------------------------------------
# Calculate the Value of a Statistical Life (following the FUND equations & 
# parameterization).
# ------------------------------------------------------------------------------

@defcomp VSL begin

    country       = Index()

    α             = Parameter(unit = "US\$2005")    # VSL scaling parameter
    ϵ             = Parameter()                     # Income elasticity of the value of a statistical life.
    y₀            = Parameter(unit = "US\$2005")    # Normalization constant.
    pc_gdp        = Parameter(index=[time, country], unit = "US\$2005/yr/person") # Country-level per capita GDP ($/person).

    gdp           = Parameter(index=[time, country], unit="billion US\$2005/yr")
    population    = Parameter(index=[time, country], unit="million")

    vsl_regions   = Parameter{Symbol}(default = :country)   # Value of a statistical life ($).
    
    vsl           = Variable(index=[time, country], unit = "US\$2005/yr")   # Value of a statistical life ($).
    
    pc_gdp_global = Variable(index=[time], unit = "US\$2005/yr/person") # Global per capita GDP ($/person).
    pc_gdp_row    = Variable(index=[time], unit = "US\$2005/yr/person") # Rest of world (row) per capita GDP ($/person).

    function run_timestep(p, v, d, t)

        ## country index 174 is currently for the USA, but this could be improved if we could sum all `except` the USA/174
        v.pc_gdp_global[t] = (sum(p.gdp[t,:])/sum(p.population[t,:])) * 1e3
        v.pc_gdp_row[t]    = ((sum(p.gdp[t,:]) - p.gdp[t,174]) / (sum(p.population[t,:]) - p.population[t,174])) * 1e3

        for c in d.country
            
            if p.vsl_regions == :country
                
                    v.vsl[t,c] = p.α * (p.pc_gdp[t,c] / p.y₀) ^ p.ϵ
                
            elseif p.vsl_regions == :us_and_row
                ## country index 174 is currently for the USA, but this iff condition could be improved by conditioning on a string "USA" in case the indices change
                if c == 174
                    v.vsl[t,c] = p.α * (p.pc_gdp[t,c] / p.y₀) ^ p.ϵ
                else
                    v.vsl[t,c] = p.α * (v.pc_gdp_row[t] / p.y₀) ^ p.ϵ
                end
                
            elseif p.vsl_regions == :global
                
                    v.vsl[t,c] = p.α * (v.pc_gdp_global[t] / p.y₀) ^ p.ϵ
                
            end
        end
    end
end
