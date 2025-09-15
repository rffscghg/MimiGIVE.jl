using Mimi

# Aggregate from countries to FUND regions using sum function

@defcomp Agriculture_RegionAggregatorSum begin

    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    input_output_mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per input region containing it's corresponding output region
    input_output_mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up

    input_region_names = Parameter{Vector{String}}(index=[ag_mapping_input_regions])
    output_region_names = Parameter{Vector{String}}(index=[ag_mapping_output_regions])

    input = Parameter(index=[time, ag_mapping_input_regions])
    output = Variable(index=[time, ag_mapping_output_regions])

    function init(p, v, d)
        idxs = indexin(p.input_output_mapping, p.output_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the Agriculture_RegionAggregatorSum's input_output_mapping Parameter must exist in the output_region_names Parameter.") : nothing
        v.input_output_mapping_int[:] = idxs
    end

    function run_timestep(p, v, d, t)
        v.output[t, :] .= 0.

        for i in d.ag_mapping_input_regions
            v.output[t, v.input_output_mapping_int[i]] += p.input[t, i]
        end
    end
end

# Version of above with no time dimension

@defcomp Agriculture_RegionAggregatorSum_NoTime begin

    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    input_output_mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per input region containing it's corresponding output region
    input_output_mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up

    input_region_names = Parameter{Vector{String}}(index=[ag_mapping_input_regions])
    output_region_names = Parameter{Vector{String}}(index=[ag_mapping_output_regions])

    input = Parameter(index=[ag_mapping_input_regions])
    output = Variable(index=[ag_mapping_output_regions])

    function init(p, v, d)
        idxs = indexin(p.input_output_mapping, p.output_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the Agriculture_RegionAggregatorSum's input_output_mapping Parameter must exist in the output_region_names Parameter.") : nothing
        v.input_output_mapping_int[:] = idxs

        # can simply fill in the data here because there is no time dimensions
        v.output[:] .= 0.
        for i in d.ag_mapping_input_regions
            v.output[v.input_output_mapping_int[i]] += p.input[i]
        end
    end

    function run_timestep(p, v, d, t)
        # blank
    end
end

# Component to disaggregate the agricultural damages from agriculture regions (FUND
# regions) to individual ISO3 countries
@defcomp AgricultureDamagesDisaggregator begin

    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    # Mapping
    mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per country containing it's corresponding region
    mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up
    fund_region_names = Parameter{Vector{String}}(index=[ag_mapping_output_regions])

    # GDP input
    gdp_country = Parameter(index=[time, ag_mapping_input_regions], unit="billion US\$2005/yr")
    gdp_fund_region = Parameter(index=[time, ag_mapping_output_regions], unit="billion US\$2005/yr")

    # Damages input
    damages_ag_fund_region = Parameter(index=[time, ag_mapping_output_regions])

    # Disaggregation
    gdp_share = Variable(index=[time, ag_mapping_input_regions]) # share of region's GDP in a given country in a given year
    damages_ag_country = Variable(index=[time, ag_mapping_input_regions])

    function init(p, v, d)
        idxs = indexin(p.mapping, p.fund_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the AgricultureDamagesDisaggregator's mapping Parameter must exist in the region_names Parameter.") : nothing
        v.mapping_int[:] = idxs
    end

    function run_timestep(p, v, d, t)

        for c in d.country
            v.gdp_share[t, c] = p.gdp_country[t, c] / p.gdp_fund_region[t, v.mapping_int[c]]
            v.damages_ag_country[t, c] = p.damages_ag_fund_region[t, v.mapping_int[c]] * v.gdp_share[t, c]
        end
    end
end
