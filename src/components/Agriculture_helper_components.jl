using Mimi

# -------------------------------------------------------------------------------------------------
# Aggregate from countries to FUND regions using sum function
# -------------------------------------------------------------------------------------------------

@defcomp Agriculture_RegionAggregatorSum begin
    
    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    input_output_mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per input region containing it's corresponding output region
    input_output_mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up

    input_region_names = Parameter{Vector{String}}(index=ag_mapping_input_regions)
    output_region_names = Parameter{Vector{String}}(index=ag_mapping_output_regions)

    input = Parameter(index=[time, ag_mapping_input_regions])
    output = Variable(index=[time, ag_mapping_output_regions])

    function init(p,v,d)
        idxs = indexin(p.input_output_mapping, p.output_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the Agriculture_RegionAggregatorSum's input_output_mapping Parameter must exist in the output_region_names Parameter.") : nothing
        v.input_output_mapping_int[:] = idxs
    end

    function run_timestep(p,v,d,t)
        v.output[t, :] .= 0.

        for i in d.ag_mapping_input_regions
            v.output[t, v.input_output_mapping_int[i]] += p.input[t,i]
        end
    end
end

# same as above but without a time dimension

@defcomp Agriculture_RegionAggregatorSum_NoTime begin
    
    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    input_output_mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per input region containing it's corresponding output region
    input_output_mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up

    input_region_names = Parameter{Vector{String}}(index=ag_mapping_input_regions)
    output_region_names = Parameter{Vector{String}}(index=ag_mapping_output_regions)

    input = Parameter(index=[ag_mapping_input_regions])
    output = Variable(index=[ag_mapping_output_regions])

    function init(p,v,d)
        idxs = indexin(p.input_output_mapping, p.output_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the Agriculture_RegionAggregatorSum's input_output_mapping Parameter must exist in the output_region_names Parameter.") : nothing
        v.input_output_mapping_int[:] = idxs

        # fill in the data because there's no time dimensions
        v.output[:] .= 0.
        for i in d.ag_mapping_input_regions
            v.output[v.input_output_mapping_int[i]] += p.input[i]
        end
    end

    function run_timestep(p,v,d,t)
        # blank
    end
end


