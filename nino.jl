#!/usr/bin/env julia
using Pkg
Pkg.add.(["HDF5","StatsBase"])
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

function find_inside_indeces(lat, lon, x0, x1, y1, y2) # Find indeces for points inside region e.g: El Nino Basin
    
    # I think this works
    x0 < 0 ? x0 += 360 : nothing
    x1 < 0 ? x1 += 360 : nothing
    y0 < 0 ? y0 += 360 : nothing
    x1 < 0 ? y1 += 360 : nothing
    
    i0 = findfirst(x -> x > x0, lon)
    i1 = findfirst(x -> x < x1, lon)
    j0 = findfirst(y -> y > y1, lat)
    j1 = findfirst(y -> y < y2, lat)

    return i0:i1, j0:j1
end

function cross_correlation(x, y)
    # Check if this looks right to you
    return (mean(x .* y) - mean(x) * mean(y)) /  std(x) / std(y)
end

function findmissing(x)
    if x < 0 # Less than 0 Kelvin
        return NaN32
    else
        return x
    end
end

function get_anomaly(;years=1948:2022, radial_period=4, scale=True)
    A = read_air_data[years[1]]
    size_A = size(A)
    close(A)
    # This is a bit silly but it should work (Allocating > 1 GB of data)
    data = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for i in 1:length(years)
        A = read_air_data(years[i])
        data[:,:,:,i] .= findmissing.(read(A["air"]))[1:365] # Discard Leap Years (:
        close(A)
    end

    out = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for i in 1:length(years)
        local_period = intersect(years, (i-radial_period):(i+radial_period))
        year_means = mean(data[:,:,:,local_period],dims=4)
        out[:,:,:,i] .-= year_means
    end
    return out 
end