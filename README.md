# MimiGIVE.jl

This package holds the scripts to run the GIVE integrated assessment model.

# 1. Preparing the Software Environment

To add the package to your current environment, run the following command at the julia package REPL:

```julia
pkg> add https://github.com/rffscghg/MimiGIVE.jl.git

```
You probably also want to install the Mimi package into your julia environment, so that you can use some of the tools from that package:

```julia
pkg> add Mimi
```

# 2. Running the Model

## 2a. Fully Coupled Model

The model uses the Mimi framework and reading the Mimi documentation first to understand the code structure is highly recommended. The basic way to access the model, run it, and explore the results is the following:

```julia
using Mimi 
using MimiGIVE

# Create the model using default specifications
m = MimiGIVE.get_model()

# Run the model
run(m)

# Explore interactive plots of all the model output.
explore(m)

# Access a specific variable's data in tabular format.
co2_emissions = m[:co2_cycle, :E_co2]

# Access a specific variable's data in a Dataframe
co2_emissions = getdataframe(m, :co2_cycle, :E_co2)
```

The `get_model` function above has the signature and options as follows:

```julia
function get_model(; Agriculture_gtap::String = "midDF",
                    socioeconomics_source::Symbol = :RFF,
                    SSP_scenario::Union{Nothing, String} = nothing,       
                    RFFSPsample::Union{Nothing, Int} = nothing,
                    Agriculture_floor_on_damages::Bool = true,
                    Agriculture_ceiling_on_benefits::Bool = false,
                    vsl::Symbol= :epa
                )
```
The relevant arguments above are described as:

**Socioeconomic**

- socioeconomics_source (default :RFF) - The options are :RFF, which uses data from 
    the RFF socioeconomic projections, or :SSP, which uses data from one of the 
    Shared Socioeconomic Pathways
    
- SSP_scenario (default to nothing) - This setting is used only if one is using 
    the SSPs as the socioeconomics_source, and the current options are "SSP119", 
    "SSP126", "SSP245", "SSP370", "SSP585", and this will be used as follows.
    See the SSPs component here: https://github.com/anthofflab/MimiSSPs.jl for more information.

    (1) Select the population and GDP trajectories for 2020 through 2300, mapping
        each RCMIP scenario to the SSP (SSP1, 2, 3, 5 respectively)
    
    (2) Choose the ar6 scenario for data from 1750 - 2019 and the RCMIP emissions 
        scenario from the MimiSSPs component to pull Leach et al. RCMIP scenario
        data for 2020 to 2300 for CO2, CH4, and N2O.

    (NOTE) that if the socioeconomics_source is :RFF this will not be consequential 
        and ssp245 will be used for the ar6 data from 1750 - 2019 and trace gases 
        from 2020 onwards, while emissions for CO2, CH4, and N2O will come from
        the MimiRFFSPs component.
    
* RFFSPsample (default to nothing, which will pull the in MimiRFFSPs) - choose
    the sample for which to run the RFF SSP. See the RFFSPs component here: 
    https://github.com/rffscghg/MimiRFFSPs.jl.

**Agriculture**

- Agriculture_gtap (default midDF) - specify the `Agriculture_gtap` input parameter as one of 
    `["AgMIP_AllDF", "AgMIP_NoNDF", "highDF", "lowDF", "midDF"]`, indicating which 
    gtap damage function the component should use. 
- Agriculture_floor_on_damages (default true) - If `Agriculture_gtap_floor_on_damages` = true, then 
    the agricultural damages (negative values of the `agcost` variable) in each 
    timestep will not be allowed to exceed 100% of the size of the  agricultural 
    sector in each region.
- Agriculture_ceiling_on_benefits (default false) - If `Agriculture_gtap_ceiling_on_benefits` = true, 
    then the agricultural benefits (positive values of the `agcost` variable) in 
    each timestep will not be allowed to exceed 100% of the size of the agricultural 
    sector in each region.
    
**Other**

- vsl (default :epa) - Specify the soruce of the value of statistical life (VSL) being used in the model. The default `:epa` uses the 2017 VSL used in U.S. EPA anlayses. Alternatively, one could use `:fund`. Both are described in the `DataExplainer` along with references to the underlying values.

## 2b. Append Offline MimiCIAM Coupling

The MimiCIAM sea level rise damages component is currently "offline" coupled to the main model due to integration barriers based on model construction and need for foresight.  To run MimiCIAM using the outputs of the GIVE model, first create and run a GIVE model as above in Section 2a., and then use the `get_ciam` and `update_ciam!` functions to obtain a CIAM model parameterized by the outputs of the GIVE model, as follows:

```julia
using Mimi 
using MimiGIVE

# Create the default GIVE model
m = MimiGIVE.get_model()

# Run the model
run(m)

# Get the default CIAM model
segment_fingerprints, m_ciam = MimiGIVE.get_ciam(m)

# Update the CIAM model with MimiGIVE specific parameters
MimiGIVE.update_ciam!(m_ciam, m, segment_fingerprints)

# Run the CIAM model
run(m_ciam)

# Access and explore the reults of CIAM using `getindex`, `getdataframe`, and
# `explore` as above

# NOTE: `explore` may be unwiedly and/or slow for CIAM, as the dimensionality 
# is large with ~12,000 coastal segments), we recommend accessing variables individually
# instead
explore(m_ciam) 
```

# 3. Running a Monte Carlo Simulation (MCS)

You can run a Monte Carlo Simulation on the GIVE model to explore the effects of parameteric uncertainty. This functionality leverages the MCS functionality of the Mimi package, and it is recommended that you review the documentation in that package for background and context.

## 3a. API

Running a basic Monte Carlo Simulation on the default model `m = MimiGIVE.get_model()` can be carried out as follows.

```julia
using Mimi 
using MimiGIVE

mcs_results = MimiGIVE.run_mcs(trials = 1000; socioeconomics_source= :RFF, output_dir = nothing, save_trials = false)
```

The built-in Monte Carlo Simulation details can be found in `src/main_mcs.jl` and the primary function has the signature as follows:

```julia
function run_mcs(;trials::Int64 = 10000, 
                    output_dir::Union{String, Nothing} = nothing, 
                    save_trials::Bool = false,
                    fair_parameter_set::Symbol = :random,
                    fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
                    rffsp_sampling::Symbol = :random, 
                    rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
                    m::Mimi.Model = get_model(), 
                    save_list::Vector = [],
                    results_in_memory::Bool = true
    )
```

This function returns the results of a Monte Carlo Simulation with the defined number of `trials`, and save data into the `output_dir` folder, optionally also saving the parameter settings for each individual trial if `save_trials` is set to `true` If no model `m` is provided we will run with the default model from `get_model()`. The rest of the arguments are described as follows: 

- `trials` (default 10,000) - number of trials to be run, used for presampling
- `output_dir` (default constructed folder name) - folder to hold results 
- `save_trials` (default false) - whether to save all random variables for all trials to trials.csv 
- `fair_parameter_set` (default :random) - :random means FAIR mcs samples will be chosen randomly from the provided sets, while :deterministic means they will be  based on the provided vector of to `fair_parameter_set_ids` keyword argument. 
- `fair_parameter_set_ids` - (default nothing) - if `fair_parameter_set` is set to :deterministic, this `n` element vector provides the fair parameter set ids that will be run, otherwise it is set to `nothing` and ignored.
- `rffsp_sampling` (default :random) - which sampling strategy to use for the RFF SPs, :random means RFF SPs will be chosen randomly, while :deterministic means they will be based on the provided vector of to `rffsp_sampling_ids` keyword argument. 
- `rffsp_sampling_ids` - (default nothing) - if `rffsp_sampling` is set to :deterministic, this `n` element vector provides the RFF SP ids that will be run, otherwise it is set to `nothing` and ignored.
- `m` (default get_model()) - the model to run the simulation for
- `save_list` (default []) - which parameters and variables to save for each trial,entered as a vector of Tuples (:component_name, :variable_name)
- `results_in_memory` (default true) - this should be turned off if you are running into memory problems, data will be streamed out to disk but not saved in memory to the mcs object

For a more in depth analysis, try exploring some saved values like temperature `T` and co2 total emissions `co2`, with 

```julia
save_list = [(:temperature, :T), (:co2_cycle, :co2)]

# default model 
mcs_results = MimiGIVE.run_mcs(trials = 100, save_list = save_list)
explore(mcs_results)

# specific model and save the trials values
m = MimiGIVE.get_model(socioeconomics_source=:SSP, SSP_scenario = "SSP585")
mcs_results = MimiGIVE.run_mcs(trials = 100, save_trials = true, m = m, save_list = save_list)
explore(mcs_results)
```

Note that if `results_in_memory` is set to `false`, you will not be able to explore your saved results from the `mcs_results` object, but instead read them in from the CSV files in the `output_dir`.

## 3b. More Details on Saving Results and Memory

**run_mcs function arguments**

- `output_dir` (default is `nothing` and then constructed under the hood) - if this is not entered, a default folder will be constructed within your `output` folder with the date, time, and trials number. Any saved data will be saved to this folder, including `trials.csv` if `save_trials = true` (a large `n` will make this too big :)) and anything in the `save_list` as described below.

- `save_list` (default is empty vector `[]`) is a vector of Tuples, each holding two Symbols, an example is below.  Note this can be any variable or parameter you see in any component of the model (all also shown in the explorer). If you are wondering about a specific parameter that is set with a random variable, take a look at the `get_mcs` function and see where a given random variable is assigned to.a component/parameter pair.  **Also feel free to inquire with the team if your curious about how to grab something specific!**
```
mcs = run_mcs(m, trials = 100, save_list = [(:temperature, :T), (:co2_cycle, :co2))
```

- `results_in_memory` (default is `true`) - if this is `true`, you will be able to access the data in the `save_list` from the returned `mcs` object with `getdataframe(mcs, :component, :variable/parameter)`, or `explore(mcs)`, but this will build up a large dataframe in memory.  If this becomes a problem, just turn this flag off as follows, and your data will **only** be streamed to files in `output_dir`.
```
mcs = run_mcs(m, trials = 100, save_list = [(:temperature, :T), (:co2_cycle, :co2), results_in_memory = false)
```

**compute_scc function arguments**

This function uses similar arguments to above, with the following differences:
- `results_in_memory` is automatically off without the option to turn it on.  This is changeable if desired.
- `output_dir` will hold **two** folders, `model1` and `model2`, which correspond to `base` and `marginal` models ie. `base` and base + pulse of gas.  This may be helpful for looking at temperature trajectories and marginal damages. Remember that the pulse units are important, take a look at these two dictionaries from `scc.jl` which may help with conversions.  **Importantly, as of now this will not include any information from CIAM damages**, you can pull anything from BRICK, but not CIAM.

```
const scc_gas_molecular_conversions = Dict(:CO2 => 12/44, # C to CO2
                                            :N2O => 28/44, # N2 to N2O,
                                            :CH4 => 1.) # CH4 to CH4
const scc_gas_pulse_size_conversions = Dict(:CO2 => 1e9, # Gt to t
                                        :N2O => 1e6, # Mt to t
                                        :CH4 => 1e6) # Mt to t
```

## 3c. Uncertain Parameters

Currently, this Monte Carlo Simulation includes the following uncertain parameters:

**Climate**
- The implementation of FAIRv1.6.2 uses the 2337 constrained parameter sets used in the AR6 (see description of details [here](https://github.com/rffscghg/MimiGIVE.jl/blob/main/docs/DataExplainer.ipynb) under the FAIR v1.6.2 heading.

- Sea Level Rise - The BRICK model varies a land water storage parameter.

**Damages**
- Socioeconomics: uncertainty using Uniform distribution across all 10,000 scenarios when RFF scenarios are enabled

- Agriculture: uncertainty using Triangular distribution across damage function parameters

- Mortality: uncertainty using resampled parameterization of damage function for Cromar et al.

- Global Damage Functions: uncertainty in Nordhaus (2017) and Howard and Sterner (2017) is derived from the parametric parameter uncertainty as stated in the corresponding publication and replication code. Kalkuhl and Wenz (2020) does not have any parameter uncertainty as it is a reduced-form approximation of their structural econometric specification. 

# 4. Calculating the SCC

We provide a user-facing API call `compute_scc` to compute the Social Cost of CO2 **in USD $\2005** for this model.  The signature of this function is as follows:

```julia
function compute_scc(m::Model=get_model(); 
                    year::Union{Int, Nothing} = nothing, 
                    last_year::Int = _model_years[end], 
                    prtp::Union{Float64,Nothing} = 0.015, 
                    eta::Union{Float64,Nothing}=1.45,
                    discount_rates=nothing,
                    certainty_equivalent=false,
                    fair_parameter_set::Symbol = :random,
                    fair_parameter_set_ids::Union{Vector{Int}, Nothing} = nothing,
                    rffsp_sampling::Symbol = :random,
                    rffsp_sampling_ids::Union{Vector{Int}, Nothing} = nothing,
                    n=0,
                    gas::Symbol = :CO2,
                    save_list::Vector = [],
                    output_dir::Union{String, Nothing} = nothing,
                    save_md::Bool = false,
                    save_cpc::Bool = false,
                    save_slr_damages::Bool = false,
                    compute_sectoral_values::Bool = false,
                    compute_domestic_values::Bool = false,
                    CIAM_foresight::Symbol = :perfect,
                    CIAM_GDPcap::Bool = false,
                    post_mcs_creation_function=nothing,
                    pulse_size::Float64=1.
    )
```

This function computes the social cost of a gas for an emissions pulse in `year` **in $2005 USD** for the provided MimiGIVE model, which can be specified as above with a particular settings. If no model is provided, the default model from MimiGIVE.get_model() is used. Furthermore, the `DamagesAggregator` component allows users to decide which damages are included in the aggregated damages used for the SCC.  The rest of the arguments are described as follows:

- `m` (default get_model()) - if no model is provided, the default model from MimiGIVE.get_model() is used. 
- 'year` (default nothing) - year of which to calculate SC (year of pulse)
- `last_year` (default 2300) - last year to include in damages summation
- `prtp` (default 0.015) and `eta` (default 1.45) - Ramsey discounting parameterization
- `discount_rates` (default nothing) - a vector of Named Tuples ie. [(prpt = 0.03., eta = 1.45), (prtp = 0.015, eta = 1.45)] - required if running n > 1
- `certainty_equivalent` (default false) - whether to compute the certainty equivalent or expected SCC
- `fair_parameter_set` (default :random) - :random means FAIR mcs samples will be chosen randomly from the provided sets, while :deterministic means they will be  based on the provided vector of to `fair_parameter_set_ids` keyword argument. 
- `fair_parameter_set_ids` - (default nothing) - if `fair_parameter_set` is set to :deterministic, this `n` element vector provides the fair parameter set ids that will be run, otherwise it is set to `nothing` and ignored.
- `rffsp_sampling` (default :random) - which sampling strategy to use for the RFF  SPs, :random means RFF SPs will be chosen randomly, while :deterministic means they will be based on the provided vector of to `rffsp_sampling_ids` keyword argument. 
- `rffsp_sampling_ids` - (default nothing) - if `rffsp_sampling` is set to :deterministic,  this `n` element vector provides the RFF SP ids that will be run, otherwise it is set to `nothing` and ignored.
- `n` (default 0) - If `n` is 0, the deterministic version will be run, otherwise, a monte carlo simulation will be run. 
- `gas` (default :CO2) - the gas for which to compute the SC, options are :CO2, :CH4, and :N2O. 
- `save_list` (default []) - which parameters and varaibles to save for each trial, entered as a vector of Tuples (:component_name, :variable_name)
- `output_dir` (default constructed folder name) - folder to hold results 
- `save_md` (default is false) - save and return the marginal damages from a monte carlo simulation
- `save_cpc` (default is false) - save and return the per capita consumption from a monte carlo simulation
- `save_slr_damages`(default is false) - save global sea level rise damages from CIAM to disk
- `compute_sectoral_values` (default is false) - compute and return sectoral values as well as total
- `compute_domestic_values` (default is false) - compute and return domestic values in addition to global
- `CIAM_foresight`(default is :perfect) - Use limited foresight (:limited) or perfect foresight (:perfect) for MimiCIAM cost calculations
- `CIAM_GDPcap` (default is false) - Limit SLR damages to country-level annual GDP
- `pulse_size` (default 1.) - This determines the size of the additional pulse of emissions. Default of `1.` implies the standard pulse size of 1Gt of C for CO2, 1Mt of CH4, and 1Mt of N2O. 

**Discount Rate Note**: In scc.jl , the rate of pure time preference `prtp` is treated as a discrete variable, such that `prtp` is applied in the discrete form `1/(1+prtp)^(t-year)`. Note that if using a `prtp` value calculated as a continuous time rate, like those in Rennert et al. (2021), one must transform this discount factor with `prtp_discrete = exp(prtp)-1)`.  Please be in touch with the developers if you need assistance or further explanation!

The social cost of a gas can be calculated either deterministically or through a Monte Carlo Simulation which incorporates parametric uncertainty and is the recommended approach. Details follow in sections 4a and 4b.

## 4a. Deterministic SCC Calculation

Running a deterministic SCC calculation uses a subset of the possible arguments to `compute_scc` to compute a single cost from default parameter specifications.

Some example use cases include:

```julia

# Compute a simple baseline case
MimiGIVE.compute_scc(year=2020)

# Compute the SCC for a different SSP/emissions scenario combination using the default sources of data (Benveniste and Leach, respectively) and a different discounting scheme parameterization
m = MimiGIVE.get_model(socioeconomics_source=:SSP, SSP_scenario="SSP585")
MimiGIVE.compute_scc(m, year=2020, prtp=0.03, eta=0.)

# Compute the partial SCC for ag:
m = MimiGIVE.get_model()
update_param!(m, :DamageAggregator, :include_ag, true)
update_param!(m, :DamageAggregator, :include_cromar_mortality, false)
update_param!(m, :DamageAggregator, :include_slr, false)

MimiGIVE.compute_scc(m, year=2020, prtp=0.03, eta=0.)
```

You can also pass `compute_scc` a vector of `NamedTuple`s to the `discount_rates` argument if you would like to compute the SCC for a few different discounting schemes.  For example:

```julia
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)]
MimiGIVE.compute_scc(m, year=2020, discount_rates = discount_rates)
```

## 4b. Monte Carlo Simulation (MCS) SCC Calculation

As described in Section 3, you can run the model using a Monte Carlo Simulation. This functionality can be used to compute a distribution of SCCs using the `n` parameter of `compute_scc`.  If `n` is 0, as is default, the deterministic SCC will be calculated.  If `n` is 2 or greater, a Monte Carlo Simulation will be run.  Note that currently for this option *you must use the `discount_rates` argument*.

The simplest case of running a Monte Carlo Simulation would look like:
```julia
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)]
result = MimiGIVE.compute_scc(year = 2020, discount_rates = discount_rates, n = 5)
```

**Optional Extra Output Specifications**

- **Marginal Damages** (only relevant for Monte Carlo Simulation): Set keyword argument `save_md` to `true` to include undiscounted marginal damages in the returned results. 
- **Net Per Capita Consumption** (only relevant for Monte Carlo Simulation): Set keyword argument `save_cpc` to `true` to include net per capita consumption in the returned results.
- **Sectorally Disaggregated Values** (only relevant for Monte Carlo Simulation): Set keyword argument `compute_sectoral_values` to `true` to compute sectorally disaggregated values. Calculations of the disaggregated sectoral SCCs will use global consumption to calculate discount factors, and thus the discount factors are consistent between the global and sectoral calculations of the SCC. To compute an isolated sectoral SCC one may run a separate simulatin with only that sector's damages turned on.
- **Within US Borders Values** - Set keyword argument `compute_domestic_values` to `true` to include SCC (and optional marginal damage) values disaggregated to the within-borders USA damages.  Calculations of the disaggregated within US borders SCC will use global consumption to calculate discount factors, and thus the discount factors are consistent between the global and within borders calculations of the SCC. 

If all four of these are set to true one would runs something like:
```julia
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)]
result = MimiGIVE.compute_scc(year = 2020, discount_rates = discount_rates, n = 5, compute_sectoral_values = true, compute_domestic_values = true, save_md = true, save_cpc = true)
```
## 4c. Returned `result` Object Structure

The object returned by `result = MimiGIVE.compute_scc(...)` is a `Dictionary` with 1-3 keys: `scc` (always), `:mds` (if `save_md` is set to `true`) and `:cpc` (if `save_cpc` is set to `true`). The structure of the values returned by these keys is as follows:
- `results[:scc]` accesses a Dictionary with keys being `NamedTuples` with elements (prtp, eta region, sector) and values which are `NamedTuples` with elements (expected_scc, se_expected_scc, and scc) as well as ce_scc and ce_sccs if certainty_equivalent=true
- `results[:mds]` accesses a Dictionary with keys being `NamedTuples` with elements (region, sector) and values which are matrices of size num trials x 281 years (2020:2300) of undiscounted marginal damages in USD $2011
- `results[:cpc]` accesses a Dictionary with keys being `NamedTuples` with elements (region, sector) and values which are matrices of size num trials x 281 years (2020:2300) of net per capita consumption in USD $2011

Below we show examples of accessing these values:

```julia
discount_rates = [(label="Ramsey", prtp=0.015, eta=1.45), (label="Constant 2%", prtp=0.02, eta=0.)]

# run the simulation with all optional outputs 
result = MimiGIVE.compute_scc(year = 2020, discount_rates = discount_rates, n = 5, compute_sectoral_values = true, compute_domestic_values = true, save_md = true, save_cpc = true)

# print out information on the calculated SCCs
for (k,v) in result[:scc]
    println("Specification: $(k.region) SCC in $(k.sector) sector, using discount rate $(k.label) specified by prtp = $(k.prtp) and eta = $(k.eta):")
    println(" --> Expected SCC = $(v.expected_scc) with standard error $(v.se_scc)")
end

# compare global and within US borders marginal damaeges
mds_global = result[:mds][(region=:globe, sector=:total)]
mds_domestic = result[:mds][(region=:domestic, sector=:total)]

# compare global and within US borders agriculture marginal damages
mds_global_ag = result[:mds][(region=:globe, sector=:agriculture)]
mds_domestic_ag = result[:mds][(region=:domestic, sector=:agriculture)]

# access net per capita consumption
net_cpc = result[:cpc][(region=:globe, sector=:total)]
```

**IMPORTANT: Please look at section 3 above for details and arguments of the Monte Carlo Simulation.**

# 5. Model Structure

Below we list the main structure of the model, for information on direct input data see docs/DataExplainer.ipynb.

## Climate

### BRICK Sea Level Rise

- Citation (1): Wong, T. E., Bakker, A. M., Ruckert, K., Applegate, P., Slangen, A., & Keller, K. (2017). BRICK v0. 2, a simple, accessible, and transparent model framework for climate and regional sea-level projections. Geoscientific Model Development, 10(7), 2741-2760.
- Citation (2): Improved Climate Modeling Reduces Extreme Social Cost of Carbon Estimates. _Submitted to Nature Climate Change_ 
- Scripts: [MimiBRICK.jl](https://github.com/raddleverse/MimiBRICK.jl)

### FAIR Climate

- Citation: Smith, C. J., Forster, P. M., Allen, M., Leach, N., Millar, R. J., Passerello, G. A., and Regayre, L. A.: FAIR v1.3: A simple emissions-based impulse response and carbon cycle model, Geosci. Model Dev., https://doi.org/10.5194/gmd-11-2273-2018, 2018.;  Millar, R. J., Nicholls, Z. R., Friedlingstein, P., and Allen, M. R.: A modified impulse-response representation of the global near-surface air temperature and atmospheric concentration response to carbon dioxide emissions, Atmos. Chem. Phys., 17, 7213-7228, https://doi.org/10.5194/acp-17-7213-2017, 2017.; specific v1.6.2 used in AR6
- Scripts: FAIR with the version noted by [IPCC-WG1](https://github.com/IPCC-WG1/Chapter-7) maps to PyPI's [FAIR version 1.6.2](https://pypi.org/project/fair/1.6.2) at the repository in [OMS-NetZero/FAIR]( https://github.com/OMS-NetZero/FAIR/tree/v1.6.2.)for the equations

### Ocean Acidification

- Citation: This package provides a simplified expression to calculate globally averaged ocean pH. Is follows Equation 7 from [Appendix F](https://www.nap.edu/read/24651/chapter/17) of "Valuing Climate Damages: Updating Estimation of the Social Cost of Carbon Dioxide": National Academies of Sciences, Engineering, and Medicine. 2017. Valuing Climate Damages: Updating Estimation of the Social Cost of Carbon Dioxide. Washington, DC: The National Academies Press.https://doi.org/10.17226/24651.
- Scripts: https://github.com/FrankErrickson/Mimi_NAS_pH.jl

## Socioeconomics (Population, GDP, and Emissions)

### RFF (CH4, CO2, and N2O)

- Citation: Rennert, K., Prest, B. C., Pizer, W. A., Newell, R. G., Anthoff, D., Kingdon, C., ... & Errickson, F. (2022). The social cost of carbon: advances in long-term probabilistic projections of population, GDP, emissions, and discount rates. Brookings Papers on Economic Activity, 2021(2), 223-305.
- Home Repository: [MimiRFFSPs.jl](https://github.com/rffscghg/MimiRFFSPs.jl)

### Shared Socioeconomic Pathways (SSPs) (CH4, CO2, SF6, and N2O)

- Citation: see MimiSSPs.jl for complete list, including Benveniste et al., 2020, Leach et al., 2021, Kikstra et al., 2021, and Riahi et al., 2017
- Home Repository: [MimiSSPs.jl](https://github.com/anthofflab/MimiSSPs.jl)

## Damages

### Sea Level Rise

- Citation: Diaz, D. B. (2016). Estimating global damages from sea level rise with the Coastal Impact and Adaptation Model (CIAM). Climatic Change, 137(1), 143-156.
- Home Repository: [MimiCIAM.jl](https://github.com/raddleverse/MimiCIAM.jl)

### Agriculture

- Citation: Moore, F.C., Baldos, U., Hertel, T. et al. New science of climate change impacts on agriculture implies higher social cost of carbon. Nat Commun 8, 1607 (2017). https://doi.org/10.1038/s41467-017-01792-x
- Home Repository: [MooreAg.jl](https://github.com/rffscghg/MooreAg.jl)

### Mortality

- Citation: Cromar, K., Howard, P., Vásquez, V. N., & Anthoff, D. (2021). Health impacts of climate change as contained in economic models estimating the social cost of carbon dioxide. GeoHealth, 5(8), e2021GH000405.

### Energy

- Citation: Clarke, L., Eom, J., Marten, E. H., Horowitz, R., Kyle, P., Link, R., ... & Zhou, Y. (2018). Effects of long-term climate change on global building energy expenditures. Energy Economics, 72, 667-677.

### Global damages functions (Kalkuhl/Wenz, Nordhaus DICE2016R2, and Howard/Sterner)

_Non-default options to use for comparisons etc._

- Citation: Kalkuhl, M., & Wenz, L. (2020). The impact of climate conditions on economic production. Evidence from a global panel of regions. Journal of Environmental Economics and Management, 103, 102360.

- Citation: Nordhaus, W. D. (2017). Revisiting the social cost of carbon. Proceedings of the National Academy of Sciences, 114(7), 1518-1523.

- Citation: Howard, P. H., & Sterner, T. (2017). Few and not so far between: a meta-analysis of climate damage estimates. Environmental and Resource Economics, 68(1), 197-225.
