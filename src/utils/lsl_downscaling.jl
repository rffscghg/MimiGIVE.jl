using Query, NetCDF, StatsBase, DataFrames, CSVFiles

# Supporting Functions to Downscale BRICK from GMSL to LSL
# Adapted from https://github.com/raddleverse/CIAM_uncertainty_propagation

"""
Retrieve BRICK fingerprints from NetCDF file
"""
function get_fingerprints(;fp_file::String = joinpath(@__DIR__, "../../data/CIAM/FINGERPRINTS_SLANGEN_Bakker.nc"))

    fplat = ncread(fp_file,"lat")
    fplon = ncread(fp_file,"lon")
    fpAIS = ncread(fp_file,"AIS")
    fpGSIC = ncread(fp_file,"GLAC")
    fpGIS = ncread(fp_file,"GIS")
    ncclose()

    return fplat,fplon,fpAIS,fpGSIC,fpGIS
end

"""
Get segment specific fingerprints for segments in segIDs_file using fingerprints in
fp_file as the baseline information. Write out these segment specific fingerprints.
"""
function get_segment_fingerprints(;fp_file::String = joinpath(@__DIR__, "../../data/CIAM/FINGERPRINTS_SLANGEN_Bakker.nc"),
                            segIDs_file::String = joinpath(@__DIR__, "../../data/CIAM/diva_segment_latlon.csv"),
                            fp_segments_file::String = joinpath(@__DIR__, "../../data/CIAM/segment_fingerprints.csv"))

    # getfingerprints from FINGERPRINTS_SLANGEN_Bakker
    (fplat,fplon,fpAIS,fpGSIC,fpGIS) = get_fingerprints(fp_file = fp_file)

    # segment data
    ciamlonlat = load(segIDs_file) |> DataFrame |> i -> sort!(i, :segments)
    ciamlonlat.longi[findall(i -> i < 0, ciamlonlat.longi)] .+= 360 # Convert Longitude to degrees East, CIAM Lat is already in (-90,90) by default

    df = DataFrame(:segments => [], :segid => [], :lon => [], :lat => [], :rgn => [],
                    :fpGIS_loc => [], 
                    :fpAIS_loc => [], 
                    :fpGSIC_loc => [], 
                    :fpTE_loc => [], 
                    :fpLWS_loc => []
    )

    for i in 1:size(ciamlonlat,1)

        lon = ciamlonlat.longi[i]
        lat = ciamlonlat.lati[i]
        segid = ciamlonlat.segid[i]
        segment = ciamlonlat.segments[i]
        rgn = ciamlonlat.rgn[i]

        # Find fingerprint degrees nearest to lat,lon
        ilat = findall(isequal(minimum(abs.(fplat.-lat))),abs.(fplat.-lat))
        ilon = findall(isequal(minimum(abs.(fplon.-lon))),abs.(fplon.-lon))

        # Take average of closest lat/lon values
        fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[ilon,ilat])))
        fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[ilon,ilat])))
        fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[ilon,ilat]))) # fixed from CIAM which had GSIC here

        fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
        fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
        fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]
        fpTE_loc = 1.0
        fpLWS_loc=1.0

        # Keep searching nearby lat/lon values if fingerprint value is NaN unless limit is hit
        inc = 1

        while isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc) && inc<5

            newlonStart = next_lon.(fplon[ilon], inc, :decrease)[1]
            newlatStart = next_lat.(fplat[ilat], inc, :decrease)[1]
            newlonEnd = next_lon.(fplon[ilon], inc, :increase)[1]
            newlatEnd = next_lat.(fplat[ilat], inc, :increase)[1]

            latInd1 = minimum(findall(isequal(minimum(abs.(fplat.-newlatStart))),abs.(fplat.-newlatStart)))
            latInd2 = maximum(findall(isequal(minimum(abs.(fplat.-newlatEnd))),abs.(fplat.-newlatEnd)))

            lonInd1 = minimum(findall(isequal(minimum(abs.(fplon.-newlonStart))),abs.(fplon.-newlonStart)))
            lonInd2 = maximum(findall(isequal(minimum(abs.(fplon.-newlonEnd))),abs.(fplon.-newlonEnd)))

            if latInd2 < latInd1
                latInds=[latInd1; 1:latInd2]
            else
                latInds=latInd1:latInd2
            end

            if lonInd2 < lonInd1
                lonInds=[lonInd1; 1:lonInd2]
            else
                lonInds = lonInd1:lonInd2
            end

            fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[lonInds,latInds])))
            fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[lonInds,latInds])))
            fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGIS[lonInds,latInds])))

            fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
            fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
            fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]

            inc = inc + 1

        end

        # If still NaN, throw an error
        if isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc)
            println("Error: no fingerprints found for ($(lon),$(lat))")
            return nothing
        end

        #append to the DataFrame
        append!(df, DataFrame(:segments => segment, :segid => segid, :lon => lon, :lat => lat, :rgn => rgn,
            :fpGIS_loc => fpGIS_loc, 
            :fpAIS_loc => fpAIS_loc, 
            :fpGSIC_loc => fpGSIC_loc, 
            :fpTE_loc => fpTE_loc, 
            :fpLWS_loc => fpLWS_loc)
        )
    end # End lonlat tuple

    df |> save(fp_segments_file)
end

"""
Downscale the data in BRICK model `m` from GMSL to LMSL using data in fp_segments_file 
as created by get_segment_fingerprints.

Output:

lsl_out: array of local sea levels, sorted in alphabetical order by segment name (time x segment)
GMSL: global mean sea levels corresponding to local sea level vector (time)
"""
function downscale_brick(m, fp_segments_file::String = joinpath(@__DIR__, "../../data/CIAM/segment_fingerprints.csv"))

    # brick data
    brick_data = DataFrame(:time => Mimi.time_labels(m),
                            :AIS => m[:global_sea_level, :slr_antartic_icesheet],
                            :GSIC => m[:global_sea_level, :slr_glaciers_small_ice_caps],
                            :GIS => m[:global_sea_level, :slr_greeland_icesheet],
                            :TE => m[:global_sea_level, :slr_thermal_expansion],
                            :LWS => m[:global_sea_level, :slr_landwater_storage],
                            :GMSL =>m[:global_sea_level, :sea_level_rise]
    )

    # segment data
    segment_fingerprints = load(fp_segments_file) |> DataFrame

    # output data
    lsl_out = zeros(size(ciamlonlat,1), size(brick_data,1)) # segments x time
    for i in 1:size(ciamlonlat,1)

       # Multiply fingerprints by BRICK ensemble members
       lsl_out[i, :] = segment_fingerprints.fpGIS_loc[i]  .* brick_data.GIS + 
                 segment_fingerprints.fpAIS_loc[i]  .* brick_data.AIS[:] + 
                 segment_fingerprints.fpGSIC_loc[i] .* brick_data.GSIC[:] +
                 segment_fingerprints.fpTE_loc[i]   .* brick_data.TE[:] + 
                 segment_fingerprints.fpLWS_loc[i]  .* brick_data.LWS[:]
                
    end # End lonlat tuple

    df = DataFrame(lsl_out, :auto) |> i -> rename!(i, Symbol.(brick_data.time)) |> DataFrame
    insertcols!(df, 1, :segid => segment_fingerprints.segid)
    insertcols!(df, 1, :segments => segment_fingerprints.segments)
    
    return df
end


##==============================================================================
## Small Helper Functions for dealing with sea level fingerprints near land

"""
    next_lat(lat::Float64, inc::Int64, direction::Symbol)
Increment latitude by `inc` in either positive direction (`direction=:increase`)
or in the negative direction (`direction=:decrease`).
Assumes latitude runs from -90 to 90 (deg N).
"""
function next_lat(lat::Float64, inc::Int64, direction::Symbol)
    if lat < -90 || lat > 90
        error("Latitude must be between -90 and 90")
    end

    if direction == :increase
        new_lat = lat + inc
        if new_lat > 90
            new_lat = new_lat - 180 #wrap around
        end

    elseif direction == :decrease
        new_lat = lat - inc
        if new_lat < -90
            new_lat = new_lat + 180
        end
    end
    return new_lat
end

"""
    next_lon(lon::Float64, inc::Int64, direction::Symbol)
Increment longitude by `inc` in either positive direction
(`direction=:increase`) or in the negative direction (`direction=:decrease`).
Assumes longitude runs from 0 to 360 (deg E).
"""
function next_lon(lon::Float64, inc::Int64, direction::Symbol)
    if lon < 0 || lon > 360
        error("Longitude must be between 0 and 360")
    end

    if direction == :increase
        new_lon = lon + inc
        if new_lon > 360
            new_lon = new_lon - 360
        end
    elseif direction == :decrease
        new_lon = lon - inc
        if new_lon < 0
            new_lon = new_lon + 360
        end
    end
    return new_lon
end
