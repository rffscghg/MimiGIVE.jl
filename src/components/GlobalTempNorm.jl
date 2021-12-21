using Mimi

# ------------------------------------------------------------------------------
# Normalize global temperature to a provided list of years - template component to be used by different components
# ------------------------------------------------------------------------------

@defcomp GlobalTempNorm begin

    global_temperature = Parameter(index=[time], unit = "degC") # Global temperature deviation (°C).

    norm_range_start = Parameter() # the first year of the range of years to normalize to
    norm_range_end = Parameter() # the last year of the range of years to normalize to

    global_temperature_norm = Variable(index=[time], unit = "degC") # Global temperature deviation normalized to the new baseline (°C).
    global_temperature_norm_range_mean = Variable(unit="degC")

    function run_timestep(p, v, d, t)

        if gettime(t) == p.norm_range_end
            t_values = TimestepValue.(collect(p.norm_range_start:1:p.norm_range_end)) # Mimi errors if you use a `:` to index with timesteps. This is a workaround for now.
            v.global_temperature_norm_range_mean = mean(p.global_temperature[t_values])
        end

        if gettime(t) >= p.norm_range_end
            v.global_temperature_norm[t] = p.global_temperature[t] - v.global_temperature_norm_range_mean
        end

    end
end
