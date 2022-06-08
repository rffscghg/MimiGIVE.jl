using MimiCIAM, Query, DataFrames, CSVFiles

# Supporting functions to support setting CIAM parameters
# Adapted from scripts in MimiCIAM.jl

"""
Process the segment-country mapping file (xsc) in CIAM by (1) Reads from CSV
and outputs list of dictionaries and arrays (2) Filters xsc file to desired
segments/regions
"""
function prep_ciam_xsc(xsc_params_path::String)

    xsc_params = load(xsc_params_path) |> DataFrame
    
    # Read in csv and convert to dictionary format
    xsc_char = Dict{Any,Any}(xsc_params.seg[i] => (xsc_params.rgn[i],xsc_params.greenland[i], xsc_params.island[i]) for i in 1:length(xsc_params.seg))

    # Create region and segment indices
    rgns = sort(unique([i[1] for i in collect(values(xsc_char))]))
    segs = string.(sort(unique(collect(keys(xsc_char)))))

    xsc_ind = Dict{Int,Tuple{Int,Int,Int}}()      # numeric seg -> (numeric rgn, greenland bool)
    xsc_segmap = Dict{Any,Any}()   # Numeric seg/rgn -> char seg/rgn
    xsc_rgnmap = Dict{Any,Any}()

    for i in 1:length(segs)
        r = xsc_char[segs[i]][1]   # Region character
        grn = xsc_char[segs[i]][2] # 0 = non-Greenland, 1 = greenland bool
        isl = xsc_char[segs[i]][3] # 0 = non-island, 1 = island bool
        r_ind = MimiCIAM.findind(r, rgns)   # Region index

        new_val = (r_ind, grn, isl)     # New tuple w/ region index instead of character

        # Build XSC Seg/rgn Maps
        r2 = rgns[r_ind]           # New region char
        s = segs[i]
        xsc_segmap[i] = s
        if !(r2 in values(xsc_rgnmap))
            xsc_rgnmap[r_ind] = r2
        end

        xsc_ind[i] = new_val
    end

    return (xsc_ind, rgns, segs, xsc_rgnmap)

end

"""
Obtain the CIAM parameters for the ciam_countries using the key in xsc_params_path
for a model with time dimension first:tstep:last and adaptation starting in `adaptation_firsts`.
"""
function get_ciam_params(;tstep::Int64, first::Int64, last::Int64, ciam_countries::Vector, xsc_params_path::String, adaptation_firsts::Array)

    # --------------------------------------------------------------------------
    # Get CIAM Default Parameters
    # Pull in main parameters and select just our countries
    ciam_params = MimiCIAM.load_ciam_params()
    for (k,v) in ciam_params 
        if "country" in names(v)
            filter!(row -> row.country in ciam_countries, ciam_params[k])
        end
    end

    # Process XSC (segment-country mapping dictionary)
    xsc_ind, rgns, segs, xsc_rgnmap  = prep_ciam_xsc(xsc_params_path)
    rgns != ciam_countries && error("The provided ciam_countries in the get_ciam_params function must match those in the provided xsc_params_path File.") : nothing

    # Process params using xsc
    MimiCIAM.parse_ciam_params!(ciam_params, rgns, segs, 0)

    # --------------------------------------------------------------------------
	# Adjust, Delete, and Add Parameters

    # --> Delete Parameters that never get used
    for p in ["s1000", "s100", "s10", "smax", "land_appr_canada", "ypc_usa", "gtapland_canada", "wbvm", "fundland_canada", "refpopdens_usa"]
        delete!(ciam_params, p)
    end

    # --> Time Related
    ciam_params["tstep"] = tstep # Length of individual time-step (years)
    ciam_params["at"] = adaptation_firsts # times that start each adaptation period
    ciam_params["ntsteps"] = length(first:tstep:last)

    # --> Metadata; not used in run
    ciam_params["rcp"] = 0
    ciam_params["percentile"] = 50
    ciam_params["ssp"] = 0

    # --> Default Settings
    ciam_params["fixed"] = true
    ciam_params["noRetreat"] = false
    ciam_params["allowMaintain"] = false
    ciam_params["popinput"] = 0
    ciam_params["discountrate"] = 0.04

    # --> IDs and Dimensions

    # Dynamically find indices corresponding to USA and CAN and manually set time steps
    # If the lengths are 0, then assume those segments are not used. Note that
    # if including Greenland, need Canada too as a reference for land appreciation

    rgn_ind_canada = [k for (k,v) in xsc_rgnmap if v=="CAN"]
    rgn_ind_canada = (length(rgn_ind_canada) > 0) ? rgn_ind_canada[1] : 0

    rgn_ind_usa = [k for (k,v) in xsc_rgnmap if v=="USA"]
    rgn_ind_usa = (length(rgn_ind_usa) > 0) ? rgn_ind_usa[1] : 0

    segID = MimiCIAM.segStr_to_segID(segs)

    ciam_params["segID"] = segID
    ciam_params["xsc"] = xsc_ind
    ciam_params["rgn_ind_canada"] = rgn_ind_canada
    ciam_params["rgn_ind_usa"] = rgn_ind_usa

    # --> Population and GDP Parameters - need to be connected to Socioeconomics

    delete!(ciam_params, "pop") # pop = Parameter(index = [time, regions])      # Population of region (million people)
    delete!(ciam_params, "ypcc") # ypcc = Parameter(index = [time, regions])     # GDP per capita per region ($2010 per capita)

    # --> Storm Damage Parameters - we adjust these to be consistent with the VSL
    # component, so remove these two parameters (see calc of vsl_ciam_country in
    # main_ciam.jl))
    delete!(ciam_params, "vslel")
    delete!(ciam_params, "vslmult")
    
    return (rgns, segs, ciam_params)
end
