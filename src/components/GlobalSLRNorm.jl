using Mimi

# ------------------------------------------------------------------------------
# Normalize global SLR to a provided list of years - template component to be used by different components
# ------------------------------------------------------------------------------

@defcomp GlobalSLRNorm begin

    global_slr = Parameter(index=[time], unit = "m")  # total sea level rise from all components (includes landwater storage for projection periods).

    norm_range_start = Parameter() # the first year of the range of years to normalize to
    norm_range_end = Parameter() # the last year of the range of years to normalize to

    global_slr_norm = Variable(index=[time], unit = "degC") # Global sea level rise deviation normalized to the new baseline (m).
    global_slr_norm_range_mean = Variable(unit="m")

    function run_timestep(p, v, d, t)

        if gettime(t) == p.norm_range_end
            t_values = TimestepValue.(collect(p.norm_range_start:1:p.norm_range_end)) # Mimi errors if you use a `:` to index with timesteps. This is a workaround for now.
            v.global_slr_norm_range_mean = mean(p.global_slr[t_values])
        end

        if gettime(t) >= p.norm_range_end
            v.global_slr_norm[t] = p.global_slr[t] - v.global_slr_norm_range_mean
        end

    end
end
