using Mimi

# Component to support summing damages across fund regions -- used to support the 
# DamageAggregator Component

@defcomp Damages_RegionAggregatorSum begin

    ag_mapping_input_regions = Index()
    ag_mapping_output_regions = Index()

    input_output_mapping = Parameter{String}(index=[ag_mapping_input_regions]) # one element per input region containing it's corresponding output region
    input_output_mapping_int = Variable{Int}(index=[ag_mapping_input_regions]) # internally computed for speed up

    input_region_names = Parameter{Vector{String}}(index=[ag_mapping_input_regions])
    output_region_names = Parameter{Vector{String}}(index=[ag_mapping_output_regions])

    damage_cromar_mortality = Parameter(index=[time, ag_mapping_input_regions], unit="US\$2005/yr")
    damage_energy = Parameter(index=[time, ag_mapping_input_regions], unit="billion US\$2005/yr")

    damage_cromar_mortality_regions = Variable(index=[time, ag_mapping_output_regions], unit="US\$2005/yr")
    damage_energy_regions = Variable(index=[time, ag_mapping_output_regions], unit="billion US\$2005/yr")

    function init(p, v, d)
        idxs = indexin(p.input_output_mapping, p.output_region_names)
        !isnothing(findfirst(i -> isnothing(i), idxs)) ? error("All provided region names in the Damages_RegionAggregatorSum's input_output_mapping Parameter must exist in the output_region_names Parameter.") : nothing
        v.input_output_mapping_int[:] = idxs
    end

    function run_timestep(p, v, d, t)
        v.damage_cromar_mortality_regions[t, :] .= 0.
        v.damage_energy_regions[t, :] .= 0.

        for i in d.ag_mapping_input_regions
            v.damage_cromar_mortality_regions[t, v.input_output_mapping_int[i]] += p.damage_cromar_mortality[t, i]
            v.damage_energy_regions[t, v.input_output_mapping_int[i]] += p.damage_energy[t, i]
        end
    end
end
