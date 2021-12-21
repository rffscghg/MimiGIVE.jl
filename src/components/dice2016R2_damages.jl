using Mimi

@defcomp dice2016R2_damage begin
    country     = Index()

    temperature = Parameter(index=[time], unit="degC")
    gdp         = Parameter(index=[time, country])
    
    a2          = Parameter(default=0.00236)

    damfrac     = Variable(index=[time])
    damages     = Variable(index=[time])

    function run_timestep(p, v, d, t)
        if p.temperature[t] < 0.
            v.damfrac[t] = 0.
        else
            v.damfrac[t] = 1 - (1/(1+(p.a2 * p.temperature[t]^2)))  # log transform to keep damages < 100%, only relevant for bad draws of the mcs. 
        end

        v.damages[t] = v.damfrac[t] * sum(p.gdp[t,:])
    end
end
