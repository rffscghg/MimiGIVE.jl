using Mimi

@defcomp hs_damage begin
    country          = Index()

    temperature      = Parameter(index=[time], unit="degC")
    gdp              = Parameter(index=[time, country])

    effects          = Parameter{Symbol}(default = :base)
    add25pct         = Parameter{Bool}(default = false)
    specification    = Parameter{Int64}(default = 7)
    
    t2_base_3        = Parameter{Float64}(default = 0.595382733860703)
    t2_prod_3        = Parameter{Float64}(default = 0.)
    t2_cat_3         = Parameter{Float64}(default = 0.259851128136597)

    t2_base_4        = Parameter{Float64}(default = 0.595382733860703)
    t2_prod_4        = Parameter{Float64}(default = 0.113324887895228)
    t2_cat_4         = Parameter{Float64}(default = 0.259851128136597)

    t2_base_7        = Parameter{Float64}(default = 0.318149737017145)
    t2_prod_7        = Parameter{Float64}(default = 0.)
    t2_cat_7         = Parameter{Float64}(default = 0.362274271711041)

    t2_base_8        = Parameter{Float64}(default = 0.318149737017145)
    t2_prod_8        = Parameter{Float64}(default = 0.398230480262918)
    t2_cat_8         = Parameter{Float64}(default = 0.362274271711041)

    t2               = Variable()
    t2_base          = Variable()
    t2_prod          = Variable()
    t2_cat           = Variable()

    damfrac          = Variable(index=[time])
    damages          = Variable(index=[time])

    function run_timestep(p, v, d, t)

        if p.specification == 3
            v.t2_base = p.t2_base_3
            v.t2_prod = p.t2_prod_3
            v.t2_cat  = p.t2_cat_3
        elseif p.specification == 4
            v.t2_base = p.t2_base_4
            v.t2_prod = p.t2_prod_4
            v.t2_cat  = p.t2_cat_4
        elseif p.specification == 7
            v.t2_base = p.t2_base_7
            v.t2_prod = p.t2_prod_7
            v.t2_cat  = p.t2_cat_7
        elseif p.specification == 8
            v.t2_base = p.t2_base_8
            v.t2_prod = p.t2_prod_8
            v.t2_cat  = p.t2_cat_8
        else
            error("Invalid effects argument of p.hs_specification")
        end
        
        ## effects options
        if p.effects == :base
            v.t2 = v.t2_base
        elseif p.effects == :productivity
            (p.specification==3 || p.specification==7 ?  error("Invalid effects argument of p.effects. This effect is not estimated in the Howard and Sterner (2017) specification p.specification") : v.t2 = v.t2_base + v.t2_prod)
        elseif p.effects == :catastrophic
            v.t2 = v.t2_base + v.t2_cat
        elseif p.effects == :total
            v.t2 = v.t2_base + v.t2_prod + v.t2_cat
        else
            error("Invalid effects argument of p.effects.")
        end

        ## 25 percent adder option
        v.t2 = p.add25pct ? v.t2*1.25 : v.t2

        ## damage function
        v.damfrac[t] = 1-(1/(1+(v.t2/100) * p.temperature[t]^2))
        v.damages[t] = v.damfrac[t] * sum(p.gdp[t,:])

    end
end
