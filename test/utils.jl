using MimiGIVE
using CSVFiles
using DataFrames
using Mimi
using Test
using Query
using Random

"""
    save_model_data(m::Model, savevars::Vector, outdir::String)

Save the data from model `m` indicated by the `savevars` vector into output 
folder `outdir`.  The `savevars` vector should hold Named Tuples (compname, varname).
"""
function save_model_data(m::Model, savevars::Vector, outdir::String)

    run(m)
    for tup in savevars
        filename = string(tup.compname, "_", tup.varname, ".csv")
        getdataframe(m, tup.compname, tup.varname) |> save(joinpath(outdir, filename))
    end
end

"""
    save_scc_data(outdir::String;
                    m::Model=MimiGIVE.get_model(),
                    year::Union{Int, Nothing} = nothing, 
                    last_year::Int = MimiGIVE._model_years[end], 
                    discount_rates = [(label="default", prtp=0.15, eta=1.45)],
                    gas::Symbol = :CO2,
                    CIAM_foresight::Symbol = :perfect,
                    CIAM_GDPcap::Bool = false,
                    pulse_size::Float64=1.
                )

Save the data from the user-specifified SCC computation using model `m` into output
folder `outdir`.
"""

function save_scc_data(outdir::String;
    m::Model=MimiGIVE.get_model(),
    year::Union{Int,Nothing}=nothing,
    last_year::Int=MimiGIVE._model_years[end],
    discount_rates=[(label="default", prtp=0.15, eta=1.45)],
    gas::Symbol=:CO2,
    CIAM_foresight::Symbol=:perfect,
    CIAM_GDPcap::Bool=false,
    pulse_size::Float64=1.
)

    df = DataFrame(:dr_label => [], :prtp => [], :eta => [], :sector => [], :scc => [])

    # global
    results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
        gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
        pulse_size=pulse_size)
    for (k, v) in results
        append!(df, DataFrame(:dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :sector => :global, :scc => v))
    end

    # check which sectors are true
    run(m)
    included_sectors = []
    for sector in [:energy, :ag, :cromar_mortality, :slr]
        if m[:DamageAggregator, Symbol(:include_, sector)]
            push!(included_sectors, sector) # add to the list
            update_param!(m, :DamageAggregator, Symbol(:include_, sector), false) # turn it off
        end
    end

    # sectoral
    for sector in included_sectors
        update_param!(m, :DamageAggregator, Symbol(:include_, sector), true)
        results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
            gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
            pulse_size=pulse_size)
        for (k, v) in results
            append!(df, DataFrame(:dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :sector => sector, :scc => v))
        end
        update_param!(m, :DamageAggregator, Symbol(:include_, sector), false)
    end

    # turn back on so m is unchanged
    for sector in included_sectors
        if m[:DamageAggregator, Symbol(:include_, sector)]
            update_param!(m, :DamageAggregator, Symbol(:include_, sector), true) # turn it on
        end
    end

    df |> save(joinpath(outdir, "SCC-$gas.csv"))
end

"""
    save_scc_mcs_data(seed::Int, outdir::String, n::Int;
                            m::Model=MimiGIVE.get_model(), 
                            year::Union{Int, Nothing} = nothing, 
                            last_year::Int = MimiGIVE._model_years[end], 
                            discount_rates=nothing,
                            certainty_equivalent=true,
                            gas::Symbol = :CO2,
                            save_list::Vector = [],
                            save_md::Bool = true,
                            save_cpc::Bool = true,
                            save_slr_damages::Bool = true,
                            compute_sectoral_values::Bool = true,
                            compute_domestic_values::Bool = true,
                            CIAM_foresight::Symbol = :perfect,
                            CIAM_GDPcap::Bool = false,
                            pulse_size::Float64 = 1.,
                        )

Save the data from the user-specifified  SCC Monte Carlo Simulation computation 
using model `m` into output folder `outdir`.
"""
function save_scc_mcs_data(seed::Int, outdir::String, n::Int;
    m::Model=MimiGIVE.get_model(),
    year::Union{Int,Nothing}=nothing,
    last_year::Int=MimiGIVE._model_years[end],
    discount_rates=nothing,
    certainty_equivalent=true,
    gas::Symbol=:CO2,
    save_list::Vector=[],
    save_md::Bool=true,
    save_cpc::Bool=true,
    save_slr_damages::Bool=true,
    compute_sectoral_values::Bool=true,
    compute_domestic_values::Bool=true,
    CIAM_foresight::Symbol=:perfect,
    CIAM_GDPcap::Bool=false,
    pulse_size::Float64=1.,
)

    Random.seed!(seed)
    results = MimiGIVE.compute_scc(m;
        n=n,
        year=year,
        last_year=last_year,
        discount_rates=discount_rates,
        certainty_equivalent=certainty_equivalent,
        fair_parameter_set=:deterministic,
        fair_parameter_set_ids=collect(1:n),
        rffsp_sampling=:deterministic,
        rffsp_sampling_ids=collect(1:n),
        gas=gas,
        save_list=save_list,
        output_dir=outdir,
        save_md=save_md,
        save_cpc=save_cpc,
        save_slr_damages=save_slr_damages,
        compute_sectoral_values=compute_sectoral_values,
        compute_domestic_values=compute_domestic_values,
        CIAM_foresight=CIAM_foresight,
        CIAM_GDPcap=CIAM_GDPcap,
        pulse_size=pulse_size
    )

    # above will save out the save list variables for model1 and model2, we now
    # need to save the scc information

    # scc
    scc_outdir = joinpath(outdir, "scc")
    mkpath(scc_outdir)

    df_expected_scc = DataFrame(:region => [], :sector => [], :dr_label => [], :prtp => [], :eta => [], :scc => [])
    df_se_expected_scc = DataFrame(:region => [], :sector => [], :dr_label => [], :prtp => [], :eta => [], :se => [])
    df_sccs = DataFrame(:region => [], :sector => [], :dr_label => [], :prtp => [], :eta => [], :scc => [], :trial => [])
    df_ce_scc = DataFrame(:region => [], :sector => [], :dr_label => [], :prtp => [], :eta => [], :scc => [])
    df_ce_sccs = DataFrame(:region => [], :sector => [], :dr_label => [], :prtp => [], :eta => [], :scc => [], :trial => [])

    for (k, v) in results[:scc]
        append!(df_expected_scc, DataFrame(:region => k.region, :sector => k.sector, :dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :scc => v.expected_scc))
        append!(df_se_expected_scc, DataFrame(:region => k.region, :sector => k.sector, :dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :se => v.se_expected_scc))
        append!(df_sccs, DataFrame(:region => k.region, :sector => k.sector, :dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :scc => v.sccs, :trial => collect(1:length(v.sccs))))
        append!(df_ce_scc, DataFrame(:region => k.region, :sector => k.sector, :dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :scc => v.ce_scc))
        append!(df_ce_sccs, DataFrame(:region => k.region, :sector => k.sector, :dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :scc => v.ce_sccs, :trial => collect(1:length(v.ce_sccs))))
    end

    df_expected_scc |> save(joinpath(scc_outdir, "expected_scc.csv"))
    df_se_expected_scc |> save(joinpath(scc_outdir, "se_expected_scc.csv"))
    df_sccs |> save(joinpath(scc_outdir, "sccs.csv"))
    df_ce_scc |> save(joinpath(scc_outdir, "ce_scc.csv"))
    df_ce_sccs |> save(joinpath(scc_outdir, "ce_sccs.csv"))

    # marginal damages
    mds_outdir = joinpath(outdir, "mds")
    mkpath(mds_outdir)

    for (k, v) in results[:mds]
        df = DataFrame(v, :auto)
        rename!(df, Symbol.(year:last_year))
        insertcols!(df, 1, :trial => 1:size(df, 1))
        df = stack(df, Not(:trial))
        df |> save(joinpath(mds_outdir, "mds-$(k.region)-$(k.sector).csv"))
    end

    # consumption per capita
    cpc_outdir = joinpath(outdir, "cpc")
    mkpath(cpc_outdir)

    region = :globe
    sector = :total

    df = DataFrame(results[:cpc][(region=region, sector=sector)], :auto)
    rename!(df, Symbol.(year:last_year))
    insertcols!(df, 1, :trial => 1:size(df, 1))
    df = stack(df, Not(:trial))
    df |> save(joinpath(cpc_outdir, "cpc-$region-$sector.csv"))

end

"""
    validate_model_data(m::Model, savevars::Vector, validationdir::String)

Validate the model `m`'s data indicated by the `savevars` vector against the 
data in the `valdiationdir`. The `savevars` vector should hold Named Tuples 
(compname, varname).
"""
function validate_model_data(m::Model, savevars::Vector, validationdir::String)

    # TOLERANCE
    rtol = 1e-9 # use relative tolerance for model data since can't assume orders of magnitude

    run(m)
    for tup in savevars

        # load validation data
        filename = string(tup.compname, "_", tup.varname, ".csv")
        validation_df = load(joinpath(validationdir, filename)) |> DataFrame

        # get the model data
        m_df = getdataframe(m, tup.compname, tup.varname)

        # test each column
        for col in names(validation_df)
            @test collect(skipmissing(validation_df[!, col])) ≈ collect(skipmissing(m_df[!, col])) rtol = rtol
        end
    end
end

"""
    validate_scc_data(validationdir::String;
                        m::Model=MimiGIVE.get_model(),
                        year::Union{Int, Nothing} = nothing, 
                        last_year::Int = MimiGIVE._model_years[end], 
                        discount_rates = [(label="default", prtp=0.15, eta=1.45)],
                        gas::Symbol = :CO2,
                        CIAM_foresight::Symbol = :perfect,
                        CIAM_GDPcap::Bool = false,
                        pulse_size::Float64=1.
                    )

Validate the the data from the user-specifified SCC computation using model `m`
against the data in the `valdiationdir`.
"""
function validate_scc_data(validationdir::String;
    m::Model=MimiGIVE.get_model(),
    year::Union{Int,Nothing}=nothing,
    last_year::Int=MimiGIVE._model_years[end],
    discount_rates=[(label="default", prtp=0.15, eta=1.45)],
    gas::Symbol=:CO2,
    CIAM_foresight::Symbol=:perfect,
    CIAM_GDPcap::Bool=false,
    pulse_size::Float64=1.
)

    # TOLERANCE
    atol = 1e-3 # for SCC dollar values

    # load validation data
    filename = "SCC-$gas.csv"
    validation_df = load(joinpath(validationdir, filename)) |> DataFrame

    # get the global model data
    results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
        gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
        pulse_size=pulse_size)

    # test each discount rate/sector combination
    for (k, v) in results
        validation_scc = validation_df |>
                         @filter(_.dr_label == k.dr_label && _.sector == "global") |>
                         DataFrame
        validation_scc = (validation_scc.scc)[1]
        @test validation_scc ≈ v atol = atol
    end

    # check which sectors are true
    run(m)
    included_sectors = []
    for sector in [:energy, :ag, :cromar_mortality, :slr]
        if m[:DamageAggregator, Symbol(:include_, sector)]
            push!(included_sectors, sector) # add to the list
            update_param!(m, :DamageAggregator, Symbol(:include_, sector), false) # turn it off
        end
    end

    for sector in included_sectors

        # get the sectoral model data
        update_param!(m, :DamageAggregator, Symbol(:include_, sector), true)
        results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
            gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
            pulse_size=pulse_size)

        # test each discount rate/sector combination
        for (k, v) in results
            validation_scc = validation_df |>
                             @filter(_.dr_label == k.dr_label && _.sector == string(sector)) |>
                             DataFrame
            validation_scc = (validation_scc.scc)[1]
            @test validation_scc ≈ v atol = atol

        end
        update_param!(m, :DamageAggregator, Symbol(:include_, sector), false)
    end

    # turn back on so m is unchanged
    for sector in included_sectors
        if m[:DamageAggregator, Symbol(:include_, sector)]
            update_param!(m, :DamageAggregator, Symbol(:include_, sector), true) # turn it on
        end
    end

end

"""
    validate_scc_mcs_data(validationdir::String, n::Int;
                            m::Model=MimiGIVE.get_model(), 
                            year::Union{Int, Nothing} = nothing, 
                            last_year::Int = MimiGIVE._model_years[end], 
                            discount_rates=nothing,
                            certainty_equivalent::Bool=true,
                            gas::Symbol = :CO2,
                            save_list::Vector = [],
                            save_md::Bool = true,
                            save_cpc::Bool = true,
                            save_slr_damages::Bool = true,
                            compute_sectoral_values::Bool = true,
                            compute_domestic_values::Bool = true,
                            CIAM_foresight::Symbol = :perfect,
                            CIAM_GDPcap::Bool = false,
                            pulse_size::Float64 = 1.
                        )

Validate the the data from the user-specifified SCC Monte Carlo Simulation computation
using model `m` against the data in the `valdiationdir`.
"""

function validate_scc_mcs_data(seed::Int, validationdir::String, n::Int;
    m::Model=MimiGIVE.get_model(),
    year::Union{Int,Nothing}=nothing,
    last_year::Int=MimiGIVE._model_years[end],
    discount_rates=nothing,
    certainty_equivalent::Bool=true,
    gas::Symbol=:CO2,
    save_list::Vector=[],
    save_md::Bool=true,
    save_cpc::Bool=true,
    save_slr_damages::Bool=true,
    compute_sectoral_values::Bool=true,
    compute_domestic_values::Bool=true,
    CIAM_foresight::Symbol=:perfect,
    CIAM_GDPcap::Bool=false,
    pulse_size::Float64=1.
)


    # TOLERANCE
    atol = 1e-3 # for SCC dollar values
    rtol = 1e-4 # use relative tolerance for non-SCC values

    # get the model data
    tmpdir = tempdir()
    Random.seed!(seed)
    results = MimiGIVE.compute_scc(m;
        n=n,
        year=year,
        last_year=last_year,
        discount_rates=discount_rates,
        certainty_equivalent=certainty_equivalent,
        fair_parameter_set=:deterministic,
        fair_parameter_set_ids=collect(1:n),
        rffsp_sampling=:deterministic,
        rffsp_sampling_ids=collect(1:n),
        gas=gas,
        save_list=save_list,
        output_dir=tmpdir,
        save_md=save_md,
        save_cpc=save_cpc,
        save_slr_damages=save_slr_damages,
        compute_sectoral_values=compute_sectoral_values,
        compute_domestic_values=compute_domestic_values,
        CIAM_foresight=CIAM_foresight,
        CIAM_GDPcap=CIAM_GDPcap,
        pulse_size=pulse_size
    )

    # save list - just compare model_1 for now, model_2 is sufficiently tested
    # by testing the scc values
    for el in save_list
        validation_df = load(joinpath(validationdir, "results", "model_1", "$(el[1])_$(el[2]).csv")) |> DataFrame
        m_df = load(joinpath(tmpdir, "results", "model_1", "$(el[1])_$(el[2]).csv")) |> DataFrame

        for col in names(validation_df) # test each column
            @test collect(skipmissing(validation_df[!, col])) ≈ collect(skipmissing(m_df[!, col])) rtol = rtol
        end
    end

    # sccs
    validation_df_expected_scc = load(joinpath(validationdir, "scc", "expected_scc.csv")) |> DataFrame
    validation_df_se_expected_scc = load(joinpath(validationdir, "scc", "se_expected_scc.csv")) |> DataFrame
    validation_df_sccs = load(joinpath(validationdir, "scc", "sccs.csv")) |> DataFrame
    validation_df_ce_scc = load(joinpath(validationdir, "scc", "ce_scc.csv")) |> DataFrame
    validation_df_ce_sccs = load(joinpath(validationdir, "scc", "ce_sccs.csv")) |> DataFrame

    for (k, v) in results[:scc]
        validation_vals = validation_df_expected_scc |>
                          @filter(_.dr_label == k.dr_label && _.region == String.(k.region) && _.sector == String.(k.sector)) |>
                          DataFrame
        validation_vals = (validation_vals.scc)[1]
        @test validation_vals ≈ v.expected_scc atol = atol

        validation_vals = validation_df_se_expected_scc |>
                          @filter(_.dr_label == k.dr_label && _.region == String.(k.region) && _.sector == String.(k.sector)) |>
                          DataFrame
        validation_vals = (validation_vals.se)[1]
        @test validation_vals ≈ v.se_expected_scc atol = atol

        validation_vals = validation_df_sccs |>
                          @filter(_.dr_label == k.dr_label && _.region == String.(k.region) && _.sector == String.(k.sector)) |>
                          DataFrame
        validation_vals = validation_vals.scc
        @test validation_vals ≈ v.sccs atol = atol

        validation_vals = validation_df_ce_scc |>
                          @filter(_.dr_label == k.dr_label && _.region == String.(k.region) && _.sector == String.(k.sector)) |>
                          DataFrame
        validation_vals = (validation_vals.scc)[1]
        @test validation_vals ≈ v.ce_scc atol = atol

        validation_vals = validation_df_ce_sccs |>
                          @filter(_.dr_label == k.dr_label && _.region == String.(k.region) && _.sector == String.(k.sector)) |>
                          DataFrame
        validation_vals = validation_vals.scc
        @test validation_vals ≈ v.ce_sccs atol = atol
    end

    # marginal damages
    for (k, v) in results[:mds]
        m_df = DataFrame(v, :auto)
        rename!(m_df, Symbol.(year:last_year))
        insertcols!(m_df, 1, :trial => 1:size(m_df, 1))
        m_df = stack(m_df, Not(:trial))

        validation_df = load(joinpath(validationdir, "mds", "mds-$(k.region)-$(k.sector).csv")) |> DataFrame
        @test validation_df[!, :value] ≈ m_df[!, :value] rtol = rtol
    end

    # consumption per capita
    region = :globe
    sector = :total

    m_df = DataFrame(results[:cpc][(region=region, sector=sector)], :auto)
    rename!(m_df, Symbol.(year:last_year))
    insertcols!(m_df, 1, :trial => 1:size(m_df, 1))
    m_df = stack(m_df, Not(:trial))

    validation_df = load(joinpath(validationdir, "cpc", "cpc-$region-$sector.csv")) |> DataFrame
    @test validation_df[!, :value] ≈ m_df[!, :value] rtol = rtol

end
