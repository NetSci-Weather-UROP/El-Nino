#!/usr/bin/env julia

using HDF5, Dates

"""
Convert NCEP dates to the actual date.
"""
function convert_time(time)
    basetime = DateTime(1800,1,1)
    time = Dates.Hour.(time) .+ basetime
end
