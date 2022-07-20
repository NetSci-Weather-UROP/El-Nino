#!/usr/bin/env julia

module OurNino

export get_data, get_anomaly, find_inside_indeces, get_period, c_i_j, in_weights

using HDF5, Dates, StatsBase, LoopVectorization

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

    j1 = findfirst(y -> y < y0, lat)-1
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
    return @turbo ((mean(x .* y) - mean(x) * mean(y)) /  std(x) / std(y))
end

function findmissing(x)
    if x < 0 # Less than 0 Kelvin
        return NaN32
    else
        return x
    end
end

function mean_air_data(A)
    size_A = size(A)
    out = Array{Float32}(undef, size_A[1], size_A[2], 365)
    @inbounds (for i in 1:size_A[1]
        for j in 1:size_A[2]
            for k in 1:365
                out[i,j,k]  = mean(findmissing.(A[i,j,(4k-3):(4k)]))
            end
        end
    end)
    return out
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
        data[:,:,:,i] .= mean_air_data(read(A["air"])) # Discard Leap Years (:
        close(A)
    end
    return data, lat, lon
end

function get_anomaly(data; years=1948:2021, radial_period=4, scale=true)
    # This is a bit silly again. We can do this on a GPU for a speedup though.
    size_A = size(data)

    out = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    @inbounds (for i in 1:length(years)
        local_period = intersect(1:length(years), (i-radial_period):(i+radial_period))
        t = Array{Task}(undef,size_A[1])
        for j in 1:size_A[1]
            t[j] = Threads.@spawn for k in 1:size_A[2]
                for l in 1:365
                    tmp = filter(x->!isnan(x),data[j,k,l,local_period])
                    out[j,k,l,i] = (data[j,k,l,i] - mean(tmp))/std(tmp)
                end
            end
        end
        wait.(t)
    end)
    return out
end

function get_period(data, y, d; days=715) #715 = 150+200+365
    size_A = size(data)
    local_data = Array{Float32}(undef, days, size_A[1], size_A[2]) 
    local_data[1:length(d:365), :, :] .= permutedims(data[:,:,d:365,y],(3,1,2))
    
    filled_count = length(d:365)
    while filled_count < days
        y+=1
        next_index = days - filled_count
        next_index > 365 ? next_index = 365 : nothing

        local_data[(filled_count+1):(filled_count+next_index), :, :] .= permutedims(data[:,:,1:next_index,y],(3,1,2))
        filled_count += next_index
    end
    return local_data
end

function c_i_j(data, is, js; lags=50:350)
    size_A = size(data)
    interior_points = length(is)*length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = Array{Float32}(undef, interior_points, exterior_points, length(lags))

    i_point_list = [(i,j) for i in is for j in js]
    e_point_list = [(i,j) for i in 1:size_A[2] for j in 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)

    t = Array{Task}(undef,interior_points*exterior_points)
    c=0
    @inbounds (for j in 1:length(e_point_list)
        for i in 1:length(i_point_list)
            c+=1
            x1 = i_point_list[i]
            x2 = e_point_list[j]
            t[c] = Threads.@spawn for d in lags .- (minimum(lags)-1)
                C[i,j, d] = cor(
                    @view(data[200:(200+364), x1[1], x1[2]]),
                    @view(data[lags[d]:(lags[d]+364), x2[1], x2[2]])
                )
            end
        end
    end)
    wait.(t)
    return C, i_point_list, e_point_list
end

function H(x)
    if x >= 0
        return typeof(x)(1)
    else
        return typeof(x)(0)
    end
end

function in_weights(C, i_point_list, e_point_list, l_lon, l_lat; lags=50:350)
    lags=200 .- lags
    size_C = size(C)
    θ = Array{Int64}(undef, size_C[1], size_C[2])
    for i in 1:size_C[1]
        for j in 1:size_C[2]
            θ[i,j] = findmax(abs,C[i,j,:])[2]
        end
    end
    
    in_C = zeros(Float32, l_lon, l_lat)

    for i in 1:size(e_point_list)[1]
        for j in 1:size(i_point_list)[1]
            
            e_point = e_point_list[i,:]

            in_C[e_point[1], e_point[2]] += (
                C[j,i,θ[j,i]]*H(lags[θ[j,i]])*(1-(i==j))
            )
        end
    end

    return in_C
end
end
