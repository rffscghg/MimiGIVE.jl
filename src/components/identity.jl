using Mimi

@defcomp IdentityComponent_co2 begin
    input_co2 = Parameter(index=[time])
    output_co2 = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output_co2[t] = p.input_co2[t]
    end
end

@defcomp IdentityComponent_n2o begin
    input_n2o = Parameter(index=[time])
    output_n2o = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output_n2o[t] = p.input_n2o[t]
    end
end

@defcomp IdentityComponent_ch4 begin
    input_ch4 = Parameter(index=[time])
    output_ch4 = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output_ch4[t] = p.input_ch4[t]
    end
end
