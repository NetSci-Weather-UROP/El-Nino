#!/usr/bin/env julia

module OurNino

export get_data,
    get_anomaly,
    get_anomaly2,
    find_inside_indeces,
    get_period,
    c_i_j,
    c_i_j2,
    c_i_j3,
    c_i_j4,
    c_i_j_full_cpu,
    in_weights

using HDF5, Dates, StatsBase, LoopVectorization

"""
Convert NCEP dates to the actual date.
"""
function convert_time(time)
    basetime = DateTime(1800, 1, 1)
    time = Dates.Hour.(time) .+ basetime
end

function read_air_data(year)
    return h5open("air.sig995/air.sig995.$year.nc")
end

function find_inside_indeces(lat, lon, x0, x1, y0, y1) # Find indeces for points inside region e.g: El Nino Basin

    # I think this works
    x0 < 0 ? x0 += 360 : nothing
    x1 < 0 ? x1 += 360 : nothing

    i0 = findfirst(x -> x > x0, lon)
    i1 = findfirst(x -> x >= x1, lon) - 1

    if x1 >= x0
        is = i0:i1
    else
        is = [1:x1; x0:length(lon)]
    end

    j1 = findfirst(y -> y <= y0, lat) - 1
    j0 = findfirst(y -> y < y1, lat)

    if y1 >= y0
        js = j0:j1
    else
        js = [1:j1; j0:length(lat)]
    end

    return is, js
end

function cross_correlation(x, y)
    # Check if this looks right to you
    return @turbo (
        (mean(x .* y) - mean(x) * mean(y)) / std(x; corrected = false) /
        std(y; corrected = false)
    )
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
    @inbounds (
        for i = 1:size_A[1]
            for j = 1:size_A[2]
                for k = 1:365
                    out[i, j, k] = mean(findmissing.(A[i, j, (4k-3):(4k)]))
                end
            end
        end
    )
    return out
end

function get_data(; years = 1948:2021)
    A = read_air_data(years[1])
    size_A = size(read(A["air"])[:, 2:(end-1), :])
    lat = read(A["lat"])[2:(end-1)]
    lon = read(A["lon"])
    close(A)
    # This is a bit silly but it should work
    data = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for i = 1:length(years)
        A = read_air_data(years[i])
        data[:, :, :, i] .= mean_air_data(read(A["air"])[:, 2:(end-1), :]) # Discard Leap Years (:
        close(A)
    end
    return data, lat, lon
end

"""
Use adjacent `radial_period` years for anomaly
"""
function get_anomaly(data; years = 1948:2021, radial_period = 4, scale = true)
    # This is a bit silly again. We can do this on a GPU for a speedup though.
    size_A = size(data)

    out = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    @inbounds (
        for i = 1:length(years)
            local_period = intersect(1:length(years), (i-radial_period):(i+radial_period))
            t = Array{Task}(undef, size_A[1])
            for j = 1:size_A[1]
                t[j] = Threads.@spawn for k = 1:size_A[2]
                    for l = 1:365
                        tmp = filter(x -> !isnan(x), data[j, k, l, local_period])
                        out[j, k, l, i] = (data[j, k, l, i] - mean(tmp)) / std(tmp)
                    end
                end
            end
            wait.(t)
        end
    )
    return out
end

"""
Use all years for anomaly as in python
"""
function get_anomaly2(data; years = 1948:2021, radial_period = 4, scale = true)
    size_A = size(data)

    out = Array{Float32}(undef, size_A[1], size_A[2], 365, length(years))
    for day = 1:365
        out[:, :, day, :] .=
            (data[:, :, day, :] .- mean(data[:, :, day, :]; dims = 3)) ./
            std(data[:, :, day, :]; dims = 3, corrected = false)
    end
    return out
end

function get_period(data, y, d; days = 715) #715 = 150+200+365
    size_A = size(data)
    local_data = Array{Float32}(undef, days, size_A[1], size_A[2])
    local_data[1:length(d:365), :, :] .= permutedims(data[:, :, d:365, y], (3, 1, 2))

    filled_count = length(d:365)
    while filled_count < days
        y += 1
        next_index = days - filled_count
        next_index > 365 ? next_index = 365 : nothing

        local_data[(filled_count+1):(filled_count+next_index), :, :] .=
            permutedims(data[:, :, 1:next_index, y], (3, 1, 2))
        filled_count += next_index
    end
    return local_data
end

"""
tau ranges between -150 and 150 here
"""
function c_i_j(data, is, js; lags = 50:350)
    size_A = size(data)
    interior_points = length(is) * length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = Array{Float32}(undef, interior_points, exterior_points, length(lags))

    i_point_list = [(i, j) for i in is for j in js]
    e_point_list = [(i, j) for i = 1:size_A[2] for j = 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)

    t = Array{Task}(undef, interior_points * exterior_points)
    c = 0
    @inbounds (
        for j = 1:length(e_point_list)
            for i = 1:length(i_point_list)
                c += 1
                x1 = i_point_list[i]
                x2 = e_point_list[j]
                t[c] = Threads.@spawn for d in lags .- (minimum(lags) - 1)
                    C[i, j, d] = cor(
                        @view(data[200:(200+364), x1[1], x1[2]]),
                        @view(data[lags[d]:(lags[d]+364), x2[1], x2[2]])
                    )
                end
            end
        end
    )
    wait.(t)
    return C, i_point_list, e_point_list
end

"""
tau ranges from 0 to 150; negative taus are interpreted to correspond to c_j_i
"""
function c_i_j2(data, is, js; lags = 150)
    size_A = size(data)
    interior_points = length(is) * length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = Array{Float32}(undef, interior_points, exterior_points, lags * 2 + 1)

    i_point_list = [(i, j) for i in is for j in js]
    e_point_list = [(i, j) for i = 1:size_A[2] for j = 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)

    t = Array{Task}(undef, interior_points * exterior_points * 2)
    c = 0
    @inbounds (
        for j = 1:length(e_point_list)
            for i = 1:length(i_point_list)
                c += 1
                x1 = i_point_list[i]
                x2 = e_point_list[j]
                t[c] = Threads.@spawn for d = 0:lags
                    C[i, j, 1+lags+d] = cor(
                        @view(data[200:(200+364), x1[1], x1[2]]),
                        @view(data[(200+d):(200+d+364), x2[1], x2[2]])
                    )
                end
                c += 1
                t[c] = Threads.@spawn for d = 1:lags
                    C[i, j, 1+lags-d] = cor(
                        @view(data[200:(200+364), x2[1], x2[2]]),
                        @view(data[(200+d):(200+d+364), x1[1], x1[2]])
                    )
                end
            end
        end
    )
    wait.(t)
    return C, i_point_list, e_point_list
end

"""
Use StatsBase' crosscor rather than implementing our own
"""
function c_i_j3(data, is, js; lags = 50:350)
    size_A = size(data)
    interior_points = length(is) * length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = Array{Float32}(undef, interior_points, exterior_points, length(lags))

    i_point_list = [(i, j) for i in is for j in js]
    e_point_list = [(i, j) for i = 1:size_A[2] for j = 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)

    t = Array{Task}(undef, interior_points * exterior_points)
    c = 0
    @inbounds (
        for j = 1:length(e_point_list)
            for i = 1:length(i_point_list)
                c += 1
                x1 = i_point_list[i]
                x2 = e_point_list[j]
                t[c] = Threads.@spawn(
                    C[i, j, :] .= crosscor(
                        @view(data[:, x1[1], x1[2]]),
                        @view(data[:, x2[1], x2[2]]),
                        lags .- 200;,
                    )
                )
            end
        end
    )
    wait.(t)
    return C, i_point_list, e_point_list
end

"""
tau ranges from 0 to 150; negative taus are interpreted to correspond to c_j_i
Computed with our cross_correlation()
"""
function c_i_j4(data, is, js; lags = 150)
    size_A = size(data)
    interior_points = length(is) * length(js)
    exterior_points = size_A[2] * size_A[3] - interior_points
    C = Array{Float32}(undef, interior_points, exterior_points, lags * 2 + 1)

    i_point_list = [(i, j) for i in is for j in js]
    e_point_list = [(i, j) for i = 1:size_A[2] for j = 1:size_A[3]]
    e_point_list = setdiff(e_point_list, i_point_list)

    t = Array{Task}(undef, interior_points * exterior_points * 2)
    c = 0
    @inbounds (
        for j = 1:length(e_point_list)
            for i = 1:length(i_point_list)
                c += 1
                x1 = i_point_list[i]
                x2 = e_point_list[j]
                t[c] = Threads.@spawn for d = 0:lags
                    C[i, j, 1+lags+d] = cross_correlation(
                        @view(data[200:(200+364), x1[1], x1[2]]),
                        @view(data[(200+d):(200+d+364), x2[1], x2[2]])
                    )
                end
                c += 1
                t[c] = Threads.@spawn for d = 1:lags
                    C[i, j, 1+lags-d] = cross_correlation(
                        @view(data[200:(200+364), x2[1], x2[2]]),
                        @view(data[(200+d):(200+d+364), x1[1], x1[2]])
                    )
                end
            end
        end
    )
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

function get_θ(C, size_C)
    θ = Array{Int64}(undef, size_C[1], size_C[2])
    for i = 1:size_C[1]
        for j = 1:size_C[2]
            θ[i, j] = findmax(abs, C[i, j, :])[2]
        end
    end
    return 0
end

function in_weights(C, i_point_list, e_point_list, l_lon, l_lat, θ; lags = 50:350)
    lags = lags .- 200

    in_C = zeros(Float32, l_lon, l_lat)

    for i = 1:size(e_point_list)[1]
        for j = 1:size(i_point_list)[1]

            e_point = e_point_list[i, :]

            in_C[e_point[1], e_point[2]] += (
                C[j, i, θ[j, i]] *
                H(lags[θ[j, i]]) *
                (i != j) *
                (abs(lags[θ[j, i]]) <= 150)
            )
        end
    end

    return in_C
end

function c_i_j_full_cpu(data; lags = (-200):200, zero_index = 200)
    size_A = size(data)
    data = permutedims(data, (2, 3, 1))
    npoints = size_A[2] * size_A[3]
    println(size_A)
    println.([size_A[2], size_A[3], size_A[2], size_A[3], length(lags)])
    C = Array{Float16}(undef, size_A[2], size_A[3], size_A[2], size_A[3], length(lags))
    means = Array{Float16}(undef, size_A[2], size_A[3], length(lags))
    mean2s = Array{Float16}(undef, size_A[2], size_A[3], length(lags))
    println("Computing Means")
    @inbounds (
        for i in axes(lags, 1)
            local_data = data[:, :, (zero_index+lags[i]):(zero_index+lags[i]+364)]
            means[:, :, i] .= mean(local_data; dims = 3)
            mean2s[:, :, i] .= mean(x -> x^2, local_data; dims = 3)
        end
    )
    println("Computing stdevs")
    stdevs = sqrt.(mean2s - (means .^ 2))
    println("Computing dots")
    @inbounds (
        for t in axes(lags, 1)
            for j in axes(C, 4)
                for i in axes(C, 2)
                    range = (zero_index):(zero_index+364)
                    range2 = (zero_index+lags[t]):(zero_index+lags[t]+364)
                    C[:, i, :, j, t] .=
                        data[:, i, range] * transpose(data[:, j, range2]) ./ 365
                end
            end
        end
    )
    println("Finishing up")
    n = findfirst(isequal(0), lags)
    @inbounds (
        for z in axes(C, 5)
            for y_alt in axes(C, 4)
                for x_alt in axes(C, 3)
                    C[:, :, x_alt, y_alt, z] .-= means[x_alt, y_alt, z] .* means[:, :, n]
                    C[:, :, x_alt, y_alt, z] ./= stdevs[x_alt, y_alt, z] .* stdevs[:, :, n]
                end
            end
        end
    )
    return C
end

end

module CudaNino

using CUDA, Main.OurNino

function cu_sums!(sums, data, lags)
    tI = threadIdx()
    bD = blockDim()
    gI = blockIdx()
    gD = gridDim()

    x_initial = tI.x
    x_step = bD.x

    y_initial = tI.y
    y_step = bD.y

    z_inital = gI.x
    z_step = gD.x

    @inbounds (
        for x = x_initial:x_step:size(data, 2)
            for y = y_initial:y_step:size(data, 3)
                for z = z_inital:z_step:size(lags, 1)
                    for i = lags[z]:(lags[z]+364)
                        sums[z, x, y] += data[200+i, x, y]
                    end
                end
            end
        end
    )

end

function cu_sum2s!(sums, data, lags)
    tI = threadIdx()
    bD = blockDim()
    gI = blockIdx()
    gD = gridDim()

    x_initial = tI.x
    x_step = bD.x

    y_initial = tI.y
    y_step = bD.y

    z_inital = gI.x
    z_step = gD.x

    @inbounds (
        for x = x_initial:x_step:size(data, 2)
            for y = y_initial:y_step:size(data, 3)
                for z = z_inital:z_step:size(lags, 1)
                    for i = lags[z]:(lags[z]+364)
                        sums[z, x, y] += data[200+i, x, y]^2
                    end
                end
            end
        end
    )
end

function cu_final_1!(C, means, data, lags, n)
    tI = threadIdx()
    bD = blockDim()
    gI = blockIdx()
    gD = gridDim()

    x_initial = tI.x
    x_step = bD.x

    y_initial = tI.y
    y_step = bD.y

    z_inital = gI.x
    z_step = gD.x

    @inbounds (
        for x = x_initial:x_step:size(data, 2)
            for y = y_initial:y_step:size(data, 3)
                for z = z_inital:z_step:size(lags, 1)
                    for x_alt in axes(data, 2)
                        for y_alt in axes(data, 3)
                            C[x, y, x_alt, y_alt, z] -=
                                means[n, x, y] * means[z, x_alt, y_alt]
                        end
                    end
                end
            end
        end
    )
end

function cu_final_2!(C, stdevs, data, lags, n)
    tI = threadIdx()
    bD = blockDim()
    gI = blockIdx()
    gD = gridDim()

    x_initial = tI.x
    x_step = bD.x

    y_initial = tI.y
    y_step = bD.y

    z_inital = gI.x
    z_step = gD.x

    @inbounds (
        for x = x_initial:x_step:size(data, 2)
            for y = y_initial:y_step:size(data, 3)
                for z = z_inital:z_step:size(lags, 1)
                    for x_alt in axes(data, 2)
                        for y_alt in axes(data, 3)
                            C[x, y, x_alt, y_alt, z] /=
                                stdevs[n, x, y] * stdevs[z, x_alt, y_alt]
                        end
                    end
                end
            end
        end
    )
end

function cu_dots!(dots, data, lags)
    tI = threadIdx()
    bD = blockDim()
    gI = blockIdx()
    gD = gridDim()

    x_initial = tI.x
    x_step = bD.x

    y_initial = tI.y
    y_step = bD.y

    z_inital = gI.x
    z_step = gD.x

    @inbounds (
        for x = x_initial:x_step:size(data, 2)
            for y = y_initial:y_step:size(data, 3)
                for z = z_inital:z_step:size(lags, 1)
                    for x_alt in axes(data, 2)
                        for y_alt in axes(data, 3)
                            for i = 0:364
                                dots[x, y, x_alt, y_alt, z] +=
                                    data[200+i, x, y] * data[200+lags[z]+i, x_alt, y_alt]
                            end
                        end
                    end
                end
            end
        end
    )
end

function c_i_j_full(data; lags = (-150):25:150)
    lags = CuArray(lags)
    size_A = size(data)
    npoints = size_A[2] * size_A[3]
    println(size_A)
    C = CuArray{Float32}(undef, size_A[2], size_A[3], size_A[2], size_A[3], length(lags))
    data = CuArray{Float32}(data)

    means = CuArray{Float32}(undef, length(lags), size_A[2], size_A[3])

    sum_kernel = @cuda launch = false cu_sums!(means, data, lags)
    sum_config = launch_configuration(sum_kernel.fun)

    println(sum_config)

    t1 = min(sum_config[:threads], size_A[2])
    t2 = min(fld(sum_config[:threads], t1), size_A[3])
    t4 = min(sum_config[:blocks], length(lags))

    sum_kernel(means, data, lags; threads = (t1, t2), blocks = t4)

    mean2s = CuArray{Float32}(undef, length(lags), size_A[2], size_A[3])

    sum2_kernel = @cuda launch = false cu_sum2s!(mean2s, data, lags)
    sum2_config = launch_configuration(sum2_kernel.fun)

    println(sum2_config)

    t1 = min(sum2_config[:threads], size_A[2])
    t2 = min(fld(sum2_config[:threads], t1), size_A[3])
    t4 = min(sum2_config[:blocks], length(lags))

    sum2_kernel(mean2s, data, lags; threads = (t1, t2), blocks = t4)

    dots_kernel = @cuda launch = false cu_dots!(C, data, lags)
    dots_config = launch_configuration(dots_kernel.fun)

    println(dots_config)

    t1 = min(dots_config[:threads], size_A[2])
    t2 = min(fld(dots_config[:threads], t1), size_A[3])
    t4 = min(dots_config[:blocks], length(lags))

    dots_kernel(C, data, lags; threads = (t1, t2), blocks = t4)

    CUDA.synchronize()

    means ./= 365
    mean2s ./= 365
    C ./= 365

    CUDA.synchronize()
    n = findfirst(isequal(0), lags)
    final_1_kernel = @cuda launch = false cu_final_1!(C, means, data, lags, n)
    final_1_config = launch_configuration(final_1_kernel.fun)

    println(final_1_config)

    t1 = min(final_1_config[:threads], size_A[2])
    t2 = min(fld(final_1_config[:threads], t1), size_A[3])
    t4 = min(final_1_config[:blocks], length(lags))

    final_1_kernel(C, means, data, lags, n; threads = (t1, t2), blocks = t4)

    stdevs = sqrt.(mean2s - (means .^ 2))

    CUDA.synchronize()

    final_2_kernel = @cuda launch = false cu_final_2!(C, means, data, lags, n)
    final_2_config = launch_configuration(final_2_kernel.fun)

    println(final_2_config)

    t1 = min(final_2_config[:threads], size_A[2])
    t2 = min(fld(final_2_config[:threads], t1), size_A[3])
    t4 = min(final_2_config[:blocks], length(lags))

    final_2_kernel(C, stdevs, data, lags, n; threads = (t1, t2), blocks = t4)

    CUDA.synchronize()

    return Array(C)
end

function full_in_weights(data, is, js, l_lon, l_lat; lags = (-150):25:150)

    in_C = zeros(Float32, l_lon, l_lat)
    @inbounds (
        for i in axes(in_C, 1)
            for j in axes(in_C, 2)
                for x in axes(in_C, 1)
                    for y in axes(in_C, 2)
                        if x in is && y in js && !(i in is && j in js)
                            θ = findmax(abs, data[x, y, i, j, :])[2]
                            #println(data[x,y,i,j,:])
                            #println(x, " ", y, " ", i, " ", j, " ", θ)
                            in_C[i, j] +=
                                (i != j) *
                                Main.OurNino.H(lags[θ]) *
                                data[x, y, i, j, θ] *
                                (abs(lags[θ]) <= 150)
                        end
                    end
                end
            end
        end
    )

    return in_C
end

end
