using Mimi

@defcomp NewSectorDamages begin

    temperature = Parameter(index=[time], unit="degC")
    gdp = Parameter(index=[time, country], unit="billion US\$2005/yr")

    a = Parameter(default=0.005)

    damfrac = Variable(index=[time, country])
    damages = Variable(index=[time, country], unit="billion US\$2005/yr")

    function run_timestep(p, v, d, t)
        if p.temperature[t] < 0.
            v.damfrac[t] = 0.
        else
            v.damfrac[t] = 1 - (1 / (1 + (p.a * p.temperature[t]^2)))
        end

        for c in d.country
            v.damages[t, c] = v.damfrac[t] * p.gdp[t, c]
        end
    end
end
