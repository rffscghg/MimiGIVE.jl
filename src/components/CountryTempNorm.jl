using Mimi

# Normalize country-level temperatures to a provided range of years.
# Template component to be used by any damage module needing country-level
# temperature normalization after pattern scaling.

@defcomp CountryTempNorm begin

    country = Index()

    country_temperature = Parameter(index=[time, country], unit="degC") # Country-level temperature (°C).

    norm_range_start = Parameter() # first year of the normalization baseline
    norm_range_end   = Parameter() # last year of the normalization baseline

    country_temperature_norm            = Variable(index=[time, country], unit="degC") # Country temperature normalized to the baseline (°C).
    country_temperature_norm_range_mean = Variable(index=[country], unit="degC")       # Per-country mean over the baseline period (°C).

    function run_timestep(p, v, d, t)

        if gettime(t) == p.norm_range_end
            t_values = TimestepValue.(collect(p.norm_range_start:1:p.norm_range_end))
            for c in d.country
                v.country_temperature_norm_range_mean[c] = mean(p.country_temperature[t_values, c])
            end
        end

        if gettime(t) >= p.norm_range_end
            for c in d.country
                v.country_temperature_norm[t, c] = p.country_temperature[t, c] - v.country_temperature_norm_range_mean[c]
            end
        end

    end
end
