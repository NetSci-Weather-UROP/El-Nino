#!/usr/bin/env julia

module OurNino

export get_data, get_anomaly, find_inside_indeces, c_i_j

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

function find_inside_indeces(lat, lon, x0, x1, y0, y1) # Find indeces for points inside region e.g: El Nino Basin
    
    # I think this works
    x0 < 0 ? x0 += 360 : nothing
    x1 < 0 ? x1 += 360 : nothing
    
    i0 = findfirst(x -> x >= x0, lon)
    i1 = findfirst(x -> x > x1, lon)-1
    
    if x1 >= x0
        is=i0:i1
    else
        is=[ 1:x1 ; x0:length(lon) ]
    end

    j1 = findfirst(y -> y < y0, lat)+1
    j0 = findfirst(y -> y <= y1, lat)

    if y1 >= y0
        js=j0:j1
    else
        js=[ 1:j1 ; j0:length(lat) ]
    end

    return is, js
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

function get_data(;years=1948:2021)
    A = read_air_data(years[1])
    size_A = size(read(A["air"]))
    lat = read(A["lat"])
    lon = read(A["lon"])
    close(A)
    # This is a bit silly but it should work
    data = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for i in 1:length(years)
        A = read_air_data(years[i])
        data[:,:,:,i] .= findmissing.(read(A["air"])[:,:,1:4:(365*4)]) # Discard Leap Years (:
        close(A)
    end
    return data, lat, lon
end

function get_anomaly(data; years=1948:2021, radial_period=4, scale=true)
    # This is a bit silly again. We can do this on a GPU for a speedup though.
    size_A = size(data)

    out = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for i in 1:length(years)
        local_period = intersect(1:length(years), (i-radial_period):(i+radial_period))
        t = Array{Task}(undef,size_A[1])
        for j in 1:size_A[1]
            t[j] = Threads.@spawn for k in 1:size_A[2]
                for l in 1:365
                    out[j,k,l,i] = data[j,k,l,i] - mean(filter(x->x!=NaN32,data[j,k,l,local_period]))
                end
            end
        end
        wait.(t)
    end
    return out
end

function get_period(data, y, d)
    local_data = Array{Float32}(undef, size_A[1], size_A[2], 365 + 200) 
    local_data[1:length(d:365)] .= data[:,:,d:365,y]
    
    filled_count = length(d:365)
    while filled_count < 365 + 200
        y+=1
        next_index = 365 + 200 - filled_count
        next_index > 365 ? nx = 365 : nothing

        local_data[(filled_count+1):(filled_count+nx)] .= data[:,:,1:nx,y]
        filled_count += next_index
    end
    return local_data
end

function c_i_j(data, is, js)
    size_A = size(data)
    interior_points = length(is)*length(js)
    exterior_points = size_A[1] * size_A[2] - interior_points
    C = Array{Array{Float32}(151)}(undef, interior_points, exterior_points)

    i_point_list = [(i,j) for i in is for j in js]
    e_point_list = [(i,j) for i in 1:size_A[1] for j in 1:size_A[2]]
    e_point_list = setdiff(e_point_list, i_point_list)

    for i in 1:length(i_point_list)
        for j in 1:length(e_point_list)
            x1 = i_point_list[i]
            x2 = e_point_list[j]
            C[i,j] = [
                cross_correlation(
                    data[x1[1], x1[2], 200:(200+364)],
                    data[x2[1], x2[2], d:(d+364)]
                ) for d in 200:-1:50
                ]
        end
    end
    return C, i_point_list, e_point_list
end
end