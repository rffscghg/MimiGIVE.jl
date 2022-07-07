using MimiGIVE
using CSVFiles
using DataFrames
using Mimi
using Test
using Query

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
                        year::Union{Int, Nothing} = nothing, 
                        last_year::Int = MimiGIVE._model_years[end], 
                        discount_rates = [(label="default", prtp=0.15, eta=1.45)],
                        gas::Symbol = :CO2,
                        CIAM_foresight::Symbol = :perfect,
                        CIAM_GDPcap::Bool = false,
                        pulse_size::Float64=1.
                    )
        
    results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
                            gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
                            pulse_size=pulse_size)
        
    df = DataFrame(:dr_label => [], :prtp => [], :eta => [], :scc => [])
    for (k,v) in results
        append!(df, DataFrame(:dr_label => k.dr_label, :prtp => k.prtp, :eta => k.eta, :scc => v))
    end
    df |> save(joinpath(outdir, "SCC-$gas.csv"))
end

"""
    save_scc_mcs_data()
"""
function save_scc_mcs_data()
    # TODO
end

"""
    validate_model_data(m::Model, savevars::Vector, validationdir::String)

Validate the model `m`'s data indicated by the `savevars` vector against the 
data in the `valdiationdir`. The `savevars` vector should hold Named Tuples 
(compname, varname).
"""
function validate_model_data(m::Model, savevars::Vector, validationdir::String)

    run(m)
    for tup in savevars
        
        # load validation data
        filename = string(tup.compname, "_", tup.varname, ".csv")
        validation_df = load(joinpath(validationdir, filename)) |> DataFrame

        # get the model data
        m_df = getdataframe(m, tup.compname, tup.varname)

        # test each column
        for col in names(validation_df)
            @test collect(skipmissing(validation_df[!, col])) ≈ collect(skipmissing(m_df[!, col])) rtol = 1e-9
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
                            year::Union{Int, Nothing} = nothing, 
                            last_year::Int = MimiGIVE._model_years[end], 
                            discount_rates = [(label="default", prtp=0.15, eta=1.45)],
                            gas::Symbol = :CO2,
                            CIAM_foresight::Symbol = :perfect,
                            CIAM_GDPcap::Bool = false,
                            pulse_size::Float64=1.
                        )
                            
    # load validation data
    filename = "SCC-$gas.csv"
    validation_df = load(joinpath(validationdir, filename)) |> DataFrame

    # get the model data
    results = MimiGIVE.compute_scc(m; year=year, last_year=last_year, discount_rates=discount_rates,
                            gas=gas, CIAM_foresight=CIAM_foresight, CIAM_GDPcap=CIAM_GDPcap,
                            pulse_size=pulse_size)
        
    # test each discount rate
    for (k,v) in results
        println(k.dr_label)
        validation_scc = validation_df |> 
            @filter(_.dr_label == k.dr_label) |> 
            DataFrame
        validation_scc = (validation_scc.scc)[1]

        @test validation_scc ≈ v atol = 1e-9
    end

end

"""
    validate_scc_mcs_data()
"""
function validate_scc_mcs_data()
    # TODO
end
