using Mimi

@defcomp AgricultureDamagesDissagregator begin

    # Mapping
    mapping = Parameter{String}(index=[country]) # one element per country containing it's corresponding region
    mapping_int = Variable{Int}(index=[country]) # internally computed for speed up
    fund_region_names = Parameter{Vector{String}}(index=[fund_region])

    damages_ag_fund_region = Parameter(index=[time, fund_region])

    # Data
    gdp_country = Parameter(index=[time, country], unit="billion US\$2005/yr")
    gdp_fund_region = Parameter(index=[time, fund_region], unit="billion US\$2005/yr")

    gdp_share   = Variable(index=[time, country])
    damages_ag_country = Variable(index=[time, country])

    function init(p,v,d)
        idxs = indexin(p.mapping, p.fund_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the AgricultureDamagesDissagregator's mapping Parameter must exist in the region_names Parameter.") : nothing
        v.mapping_int[:] = idxs
    end

    function run_timestep(p,v,d,t)

        for c in d.country
            v.gdp_share[t,c] = p.gdp_country[t, c] / p.gdp_fund_region[t, mapping_int[c]]
            v.damages_ag_country = p.damges_ag_fund_region[t, mapping_int[c]] * v.gdp_share[t,c]
        end

    end
end
