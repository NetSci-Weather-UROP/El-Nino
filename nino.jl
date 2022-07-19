#!/usr/bin/env julia

module OurNino

export get_data, get_anomaly, find_inside_indeces, get_period, c_i_j

using HDF5, Dates, StatsBase, LoopVectorization, CUDA

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
                out[i,j,k] = mean(filter(x->x>0,A[i,j,(4k-3):(4k)]))
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
    return CuArray(data), lat, lon
end

function remove_anomaly!(out, data, s_data, l_years, radial_period)
    ti=threadIdx()
    x=ti.x
    y=ti.y
    T=ti.z

    z=blockIdx().x
    zm = blockDim().x

    tm = gridDim()
    xm = tm.x
    ym = tm.y
    Tm = tm.z
    
    @inbounds(for l in z:zm:l_years
        ll=max(1,l-radial_period)
        ul=min(l_years,l+radial_period)
        n=ul-ll+1
            for k in T:Tm:365
                for j in y:ym:s_data[2]
                    for i in x:xm:s_data[1]
                    out[i,j,k,l] = data[i,j,k,l] - sum(@view(data[x,y,z,ll:ul]))/n
                end
            end
        end
    end)
    return
end

function get_anomaly(data; years=1948:2021, radial_period=4, scale=true)
    size_A = size(data)

    out = CuArray{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    
    kernel = @cuda launch=false remove_anomaly!(out, data, size_A, length(years), radial_period)
    config = launch_configuration(kernel.fun)

    t1 = min(size_A[1],config[2])
    t2 = min(fld(config[2],t1),size_A[2])
    t3 = min(fld(config[2],t1*t2),length(years))
    t4 = min(365,config[1])

    println(config)
    println(size_A)
    println(t1, " ", t2, " ", t3, " ", t4)

    kernel(out, data, size_A, length(years), radial_period; threads=(t1,t2,t3), blocks=t4)
    return out
end

function get_period(data, y, d; days=715) #715 = 150+200+365
    size_A = size(data)
    local_data = CuArray{Float32}(undef, days, size_A[1], size_A[2]) 
    local_data[1:length(d:365), :, :] .= permutedims(@view(data[:,:,d:365,y]),(3,1,2))
    
    filled_count = length(d:365)
    while filled_count < days
        y+=1
        next_index = days - filled_count
        next_index > 365 ? next_index = 365 : nothing

        local_data[(filled_count+1):(filled_count+next_index), :, :] .= permutedims(@view(data[:,:,1:next_index,y]),(3,1,2))
        filled_count += next_index
    end
    return local_data
end

function xsq(x)
    return x^2
end

function c_i_j!(data, C, i_point_list, e_point_list, interior_points, exterior_points, lags, l_lags)
    ti=threadIdx()
    x=ti.x
    y=ti.y

    z=blockIdx().x
    zm = blockDim().x

    tm = gridDim()
    xm = tm.x
    ym = tm.y

    @inbounds(for j in y:ym:exterior_points
        x2 = @view(e_point_list[j, :])
        for i in x:xm:interior_points
            x1 = @view(i_point_list[i, :])
            for d in z:zm:l_lags
                dotmean = 0f0
                for day in 0:364
                    dotmean += data[(200+day), x1[1], x1[2]] * data[lags[d]+day, x2[1], x2[2]]
                end
                dotmean /= 365
                mean1 = sum(@view(data[200:(200+364), x1[1], x1[2]]))/365
                mean2 = sum(@view(data[lags[d]:(lags[d]+364), x2[1], x2[2]]))/365
                meansq1 = sum(xsq,@view(data[200:(200+364), x1[1], x1[2]]))/365
                meansq2 = 0f0
                for day in 0:364
                    meansq2 += data[lags[d]+day,x2[1],x2[2]]^2
                end
                meansq2 /= 365 
                std1 = sqrt(meansq1 - mean1^2)
                std2 = sqrt(meansq2 - mean2^2)
                
                C[i,j,d] = (dotmean - mean1 * mean2) / std1 / std2
            end
        end
    end)
    return
end

function c_i_j(data, is, js; lags=50:350)
    size_A = size(data)
    interior_points = length(is)*length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = CuArray{Float32}(undef, interior_points, exterior_points, length(lags))

    i_point_list = [(i,j) for i in is for j in js]
    e_point_list = [(i,j) for i in 1:size_A[2] for j in 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)
    i_point_list = CuArray(hcat(first.(i_point_list),last.(i_point_list)))
	e_point_list = CuArray(hcat(first.(e_point_list),last.(e_point_list)))

    kernel = @cuda launch=false c_i_j!(data, C, i_point_list, e_point_list, interior_points, exterior_points, lags, length(lags))
    config = launch_configuration(kernel.fun)

    t1 = min(interior_points,config[2])
    t2 = min(fld(config[2],t1),exterior_points)
    t4 = min(length(lags),config[1])

    println(config)
    println(t1, " ", t2, " ", t4)
    
    kernel(data, C, i_point_list, e_point_list, interior_points, exterior_points, lags, length(lags); threads=(t1,t2), blocks=t4)
    return C, i_point_list, e_point_list
end
end
