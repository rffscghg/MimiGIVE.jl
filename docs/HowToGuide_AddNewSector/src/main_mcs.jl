using Mimi
using MimiGIVE
using Distributions
using Dates

# used for Option 1 in Advanced Topic: Uncertainty and Intermediate Outputs
function update_mcs!(mcs)

    # add new sector uncertainty
    rv_name = :rv_new_sector_a
    Mimi.add_RV!(mcs, rv_name, Normal(0.005, 0.005 / 2)) # add random variable
    Mimi.add_transform!(mcs, :NewSectorDamages, :a, :(=), rv_name) # connect random variable to parameter

end

# used for Option 2 in Advanced Topic: Uncertainty and Intermediate Outputs
function get_modified_mcs(trials; args...)
    mcs = MimiGIVE.get_mcs(trials; args...) # get the original MCS

    # add new sector uncertainty
    rv_name = :rv_new_sector_a
    Mimi.add_RV!(mcs, rv_name, Normal(0.005, 0.005 / 2)) # add random variable
    Mimi.add_transform!(mcs, :NewSectorDamages, :a, :(=), rv_name) # connect random variable to parameter

    return mcs
end

function run_modified_mcs(; trials::Int64=10000,
    output_dir::Union{String,Nothing}=nothing,
    save_trials::Bool=false,
    fair_parameter_set::Symbol=:random,
    fair_parameter_set_ids::Union{Vector{Int},Nothing}=nothing,
    rffsp_sampling::Symbol=:random,
    rffsp_sampling_ids::Union{Vector{Int},Nothing}=nothing,
    m::Mimi.Model=get_modified_model(), # <-- using a different default model
    save_list::Vector=[],
    results_in_memory::Bool=true,
)

    m = deepcopy(m) # in the case that an `m` was provided, be careful that we don't modify the original

    trials < 2 && error("Must run `run_mcs` function with a `trials` argument greater than 1 due to a Mimi specification about SampleStores.  TO BE FIXED SOON!")

    # Set up output directories
    output_dir = output_dir === nothing ? joinpath(@__DIR__, "../output/mcs/", "MCS $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$trials") : output_dir
    isdir("$output_dir/results") || mkpath("$output_dir/results")

    trials_output_filename = save_trials ? joinpath("$output_dir/trials.csv") : nothing

    socioeconomics_module = MimiGIVE._get_module_name(m, :Socioeconomic)
    if socioeconomics_module == :MimiSSPs
        socioeconomics_source = :SSP
    elseif socioeconomics_module == :MimiRFFSPs
        socioeconomics_source = :RFF
    end

    # Get an instance of the mcs
    mcs = get_modified_mcs(trials;  # <-- using a different function to obtain the mcs
        socioeconomics_source=socioeconomics_source,
        mcs_years=Mimi.time_labels(m),
        fair_parameter_set=fair_parameter_set,
        fair_parameter_set_ids=fair_parameter_set_ids,
        rffsp_sampling=rffsp_sampling,
        rffsp_sampling_ids=rffsp_sampling_ids,
        save_list=save_list,
    )

    # run monte carlo trials
    results = run(mcs,
        m,
        trials;
        trials_output_filename=trials_output_filename,
        results_output_dir="$output_dir/results",
        results_in_memory=results_in_memory
    )

    return results
end
