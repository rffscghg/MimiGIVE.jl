using Mimi

# ------------------------------------------------------------------------------
# Accumulate the Ocean Heat Content from BRICK Over Time
# ------------------------------------------------------------------------------

@defcomp OceanHeatAccumulator begin

    del_ohc         = Parameter(index=[time], unit="J") # year over year Ocean heat content anomaly
    del_ohc_accum   = Variable(index=[time], unit="1e22 J") # accumulated Ocean heat content anomaly

    function run_timestep(p, v, d, t)

        # The BRICK TE component multiplies ocean heat by 1e22 because it's assuming 
        # the SNEASY units. FAIR ocean heat is already in units of 10^22, so this
        # divides by 1e22 just so it can be re-scaled again in the BRICK TE component.
        if is_first(t)
            v.del_ohc_accum[t] = 0. # FAIR won't provide del_ohc for first period so leave at 0.
        else
            v.del_ohc_accum[t] = v.del_ohc_accum[t-1] + (p.del_ohc[t] ./ 1e22)
        end
    end
end
