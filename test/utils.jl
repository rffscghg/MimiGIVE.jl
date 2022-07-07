using MimiGIVE
using CSVFiles
using DataFrames
using Mimi
using Test

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
    save_scc_data(m::Model, outdir::String)

Save the data from an scc computation from model `m` indicated into output folder 
`outdir`.  
"""
function save_scc_data(m::Model, outdir::String)
    # TODO
end

"""
    save_scc_mcs_data()
"""
function save_scc_mcs_data()
    # TODO
end

"""
Validate the model `m`'s data indicated by the `savevars` vector against the 
data in the `valdiationdir`. The `savevars` vector should hold Named Tuples 
(compname, varname).

    validate_model_data(m::Model, savevars::Vector, validationdir::String)
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
    validate_scc_data(m::Model, validationdir::String)  

Validate the model `m`'s scc computation data indicated against the 
data in the `valdiationdir`.
"""
function validate_scc_data(m::Model, validationdir::String)
    # TODO
end

"""
    validate_scc_mcs_data()
"""
function validate_scc_mcs_data()
    # TODO
end
