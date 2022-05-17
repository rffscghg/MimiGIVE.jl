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

    vsl           = Variable(index=[time, country], unit = "US\$2005/yr")   # Value of a statistical life ($).
    
    function run_timestep(p, v, d, t)
        for c in d.country
            v.vsl[t,c] = p.α * (p.pc_gdp[t,c] / p.y₀) ^ p.ϵ
        end
    end
end
