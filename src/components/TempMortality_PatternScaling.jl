using Mimi

# ------------------------------------------------------------------------------
# Calculate population-weighted, country-level temperatures which are provided relative to 2001-2020 mean
# ------------------------------------------------------------------------------

@defcomp TempMortality_PatternScaling begin

    country = Index()
    cmip6_gcms = Index()

    gcm_id             = Parameter{Int64}(default = 5) # the pattern gcm id to use, default to CNRM-CM6-1
    pattern            = Parameter(index=[country, cmip6_gcms]) # Population or GDP-weighted coefficients that scale global to local temperatures (Δ country °C / Δ global °C)
    global_temperature = Parameter(index=[time], unit = "degC") # Global average surface temperature anomaly (°C). - should be normalized to 2001 - 2020

    local_temperature = Variable(index=[time,country], unit = "degC") # Country-level temperatures derived from the pattern (°C).

    function run_timestep(p, v, d, t)
        for c in d.country
            # Calculate population-weighted country-level temperatures.
            v.local_temperature[t,c] = p.pattern[c, p.gcm_id] * p.global_temperature[t]
        end
    end
end
