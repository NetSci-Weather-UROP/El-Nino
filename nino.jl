#!/usr/bin/env julia

using HDF5, Dates, StatsBase

"""
Convert NCEP dates to the actual date.
"""
function convert_time(time)
    basetime = DateTime(1800,1,1)
    time = Dates.Hour.(time) .+ basetime
end

function read_air_data(year)
    return h5open("air.sig995/air.sig995.$year.nc")
end

struct pos_data # Positional data
    data::Matrix
    lat::Vector
    lon::Vector
end

function find_inside(data::pos_data, x0, x1, y1, y2) # Find points inside region e.g: El Nino Basin
    
    # I think this works
    x0 < 0 ? x0 += 360 : nothing
    x1 < 0 ? x1 += 360 : nothing
    y0 < 0 ? y0 += 360 : nothing
    x1 < 0 ? y1 += 360 : nothing
    
    i0 = findfirst(
        x -> x > x0,
        data.lon
        )
    
    i1 = findfirst(
        x -> x < x1,
        data.lon
    )

    j0 = findfirst(
        y -> y > y1,
        data.lat
    )
    
    j1 = findfirst(
        y -> y < y2,
        data.lat
    )

    return pos_data(
        data.data[i0:i1, j0:j1],
        data.lat[i0:i1],
        data.lon[j0:j1]
    )
end