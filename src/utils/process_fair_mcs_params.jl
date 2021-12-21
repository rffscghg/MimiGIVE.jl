using JSON, CSVFiles, DataFrames

# Process all FAIR Monte Carlo Simulation constrained parameter sets from JSON
# to CSV files to use directly by functions in main_mcs.jl

# Pairings of parameters to csv files names, for reference, also produced by
# `get_fair_mcs_params_map` in `main_mcs.jl`

# Dict(
#     :β_CO                          => "b_aero_CO",
#     :scale_CH₄                     => "scale_CH4",
#     :F_solar                       => "F_solar",
#     :Ψ_CH₄                         => "b_tro3_CH4",
#     :scale_N₂O                     => "scale_N2O",
#     :CO₂_pi                        => "C_pi",
#     :deep_ocean_efficacy           => "deep_ocean_efficacy",
#     :scale_bcsnow                  => "scale_bcsnow",
#     :scale_aerosol_direct_OC       => "scale_aerosol_direct_OC",
#     :b_SOx                         => "ghan_params_SOx",
#     :feedback                      => "ozone_feedback",
#     :scale_O₃                      => "scale_O3",
#     :b_POM                         => "ghan_params_b_POM",
#     :r0_co2                        => "r0",
#     :β_NH3                         => "b_aero_NH3",
#     :lambda_global                 => "lambda_global",
#     :scale_landuse                 => "scale_landuse",
#     :scale_volcanic                => "scale_volcanic",
#     :scale_aerosol_direct_SOx      => "scale_aerosol_direct_SOx",
#     :β_NOx                         => "b_aero_NOx",
#     :Ψ_N₂O                         => "b_tro3_N2O",
#     :ocean_heat_capacity           => "ocean_heat_capacity",
#     :β_OC                          => "b_aero_OC",
#     :scale_solar                   => "scale_solar",
#     :rC_co2                        => "rc",
#     :scale_aerosol_direct_BC       => "scale_aerosol_direct_BC",
#     :scale_CH₄_H₂O                 => "scale_CH4_H2O",
#     :scale_aerosol_indirect        => "scale_aerosol_indirect",
#     :scale_ods                     => "scale_ods",
#     :Ψ_CO                          => "b_tro3_CO",
#     :scale_aerosol_direct_NOx_NH3  => "scale_aerosol_direct_NOx_NH3",
#     :scale_other_ghg               => "scale_other_ghg",
#     :Ψ_NMVOC                       => "b_tro3_NMVOC",
#     :F2x                           => "F2x",
#     :β_SOx                         => "b_aero_SOx",
#     :β_NMVOC                       => "b_aero_NMVOC",
#     :rT_co2                        => "rt",
#     :β_BC                          => "b_aero_BC",
#     :scale_CO₂                     => "scale_CO2",
#     :Ψ_ODS                         => "b_tro3_ODS",
#     :scale_aerosol_direct_CO_NMVOC => "scale_aerosol_direct_CO_NMVOC",
#     :Ψ_NOx                         => "b_tro3_NOx",
#     :ocean_heat_exchange           => "ocean_heat_exchange",
#     :ϕ                             => "ghan_params_Pi"
# )

n = 2237 # total number of available samples
fair_params = JSON.parsefile(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair-1.6.2-wg3-params.json"));

# Names of minor greenhouse gases and ozone-depleting substances (used or indexing).
other_ghg_names = ["CF4", "C2F6", "C6F14", "HFC23", "HFC32", "HFC43_10", "HFC125", "HFC134a", "HFC143a", "HFC227ea", "HFC245fa", "SF6"]
ods_names       = ["CFC_11", "CFC_12", "CFC_113", "CFC_114", "CFC_115", "CARB_TET", "MCF", "HCFC_22", "HCFC_141B", "HCFC_142B", "HALON1211", "HALON1202", "HALON1301", "HALON2402", "CH3BR", "CH3CL"]

# Carbon cycle

for p in ["r0", "rt", "rc"]
    DataFrame(p => [fair_params[i][p] for i in 1:n]) |> 
        save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))
end

# Forcing from a doubling of CO₂.

for p in ["F2x"]
    DataFrame(p => [fair_params[i][p] for i in 1:n]) |> 
        save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))
end

# Ozone radiative forcing feedback.

for p in ["ozone_feedback"]
    DataFrame(p => [fair_params[i][p] for i in 1:n]) |> 
        save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))
end

# Solar radiative forcing. 

for p in ["C_pi"] # Choose first element of vector of 31 elements
    DataFrame(p => [fair_params[i][p][1] for i in 1:n]) |>
        save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))
end

# Pre-industrial CO₂ concentration (other concentrations fixed across samples). 

for p in ["F_solar"]
    arr = [fair_params[i][p] for i in 1:n]
    arr = reduce(hcat, arr)' # 361 years per sample - flatten out from vector of vectors to a matrix
    df = DataFrame(arr, :auto) |>
    i -> rename!(i, Symbol.(1750:2110)) # TODO are these the correct years?
    df |> save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))
end

# Temperature component

for p in ["ocean_heat_exchange", "deep_ocean_efficacy", "lambda_global", "ocean_heat_capacity"]

    if p == "ocean_heat_capacity"
        arr = [fair_params[i][p] for i in 1:n]
        arr = reduce(hcat, arr)' # 2 members (deep and mixed) per sample - flatten out from vector of vectors to a matrix
        df = DataFrame(arr, :auto)
        rename!(df, ["1", "2"])
    else
        df = DataFrame(p => [fair_params[i][p] for i in 1:n])
    end
    df |> save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_$p.csv"))

end

# "ghan_params" for aerosol indirect forcing effect.

DataFrame(:ghan_params_Pi => [fair_params[i]["ghan_params"][1] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_ghan_params_Pi.csv"))

DataFrame(:ghan_params_SOx => [fair_params[i]["ghan_params"][2] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_ghan_params_SOx.csv"))

DataFrame(:ghan_params_b_POM => [fair_params[i]["ghan_params"][3] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_ghan_params_b_POM.csv"))

# Radiative forcing scaling terms (based on ordering of forcing agents in Python code). - select from a vector of 45 elements

# TODO: :scale_contrails !!! Default FAIR has contrail forcing switched off. But they sample a scaling term. Not including for now.

DataFrame(:scale_CO2 => [fair_params[i]["scale"][1] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_CO2.csv"))

DataFrame(:scale_CH4 => [fair_params[i]["scale"][2] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_CH4.csv"))

DataFrame(:scale_N2O => [fair_params[i]["scale"][3] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_N2O.csv"))

DataFrame(:scale_O3 => [fair_params[i]["scale"][32] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_O3.csv"))

DataFrame(:scale_CH4_H2O => [fair_params[i]["scale"][34] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_CH4_H2O.csv"))

DataFrame(:scale_aerosol_direct_SOx => [fair_params[i]["scale"][36] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_direct_SOx.csv"))

DataFrame(:scale_aerosol_direct_CO_NMVOC => [fair_params[i]["scale"][37] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_direct_CO_NMVOC.csv"))

DataFrame(:scale_aerosol_direct_NOx_NH3 => [fair_params[i]["scale"][38] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_direct_NOx_NH3.csv"))

DataFrame(:scale_aerosol_direct_BC => [fair_params[i]["scale"][39] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_direct_BC.csv"))

DataFrame(:scale_aerosol_direct_OC => [fair_params[i]["scale"][40] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_direct_OC.csv"))

DataFrame(:scale_aerosol_indirect => [fair_params[i]["scale"][41] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_aerosol_indirect.csv"))

DataFrame(:scale_bcsnow => [fair_params[i]["scale"][42] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_bcsnow.csv"))

DataFrame(:scale_landuse => [fair_params[i]["scale"][43] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_landuse.csv"))

DataFrame(:scale_volcanic => [fair_params[i]["scale"][44] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_volcanic.csv"))

DataFrame(:scale_solar => [fair_params[i]["scale"][45] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_solar.csv"))

scale_other_ghg = [fair_params[i]["scale"][4:15] for i in 1:n]
scale_other_ghg = reduce(hcat, scale_other_ghg)'
scale_other_ghg = DataFrame(scale_other_ghg, :auto) |>
    i -> rename!(i, other_ghg_names) |>
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_other_ghg.csv"))

scale_ods = [fair_params[i]["scale"][16:31] for i in 1:n]
scale_ods = reduce(hcat, scale_ods)'
scale_ods = DataFrame(scale_ods, :auto) |>
    i -> rename!(i, ods_names) |>
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_scale_ods.csv"))

# Ozone radiative forcing - select from a vector of 6 elements

DataFrame(:b_tro3_CH4 => [fair_params[i]["b_tro3"][1] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_CH4.csv"))

DataFrame(:b_tro3_N2O => [fair_params[i]["b_tro3"][2] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_N2O.csv"))

DataFrame(:b_tro3_ODS => [fair_params[i]["b_tro3"][3] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_ODS.csv"))

DataFrame(:b_tro3_CO => [fair_params[i]["b_tro3"][4] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_CO.csv"))

DataFrame(:b_tro3_NMVOC => [fair_params[i]["b_tro3"][5] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_NMVOC.csv"))

DataFrame(:b_tro3_NOx => [fair_params[i]["b_tro3"][6] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_tro3_NOx.csv"))

# Aerosol direct forcing.

DataFrame(:b_aero_SOx => [fair_params[i]["b_aero"][1] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_SOx.csv"))

DataFrame(:b_aero_CO => [fair_params[i]["b_aero"][2] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_CO.csv"))

DataFrame(:b_aero_NMVOC => [fair_params[i]["b_aero"][3] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_NMVOC.csv"))

DataFrame(:b_aero_NOx => [fair_params[i]["b_aero"][4] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_NOx.csv"))

DataFrame(:b_aero_BC => [fair_params[i]["b_aero"][5] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_BC.csv"))

DataFrame(:b_aero_OC => [fair_params[i]["b_aero"][6] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_OC.csv"))

DataFrame(:b_aero_NH3 => [fair_params[i]["b_aero"][7] for i in 1:n]) |> 
    save(joinpath(@__DIR__, "..", "..", "data", "FAIR_mcs", "fair_mcs_params_b_aero_NH3.csv"))
