#!/usr/bin/env julia

using HDF5, CairoMakie, GeoMakie

function H(x)
    if x >= 0
        return typeof(x)(1)
    else
        return typeof(x)(0)
    end
end

function main()
    lags = 200 .- (50:350)

    generated_data = h5open("generated_data.h5", "r")

    C = read(generated_data["C"])
    size_C = size(C)

    lat = read(generated_data["lat"])
    lon = read(generated_data["lon"])

    e_point_list = read(generated_data["e_point_list"])
    i_point_list = read(generated_data["i_point_list"])

    θ = Array{Int64}(undef, size_C[1], size_C[2])
    for i = 1:size_C[1]
        for j = 1:size_C[2]
            θ[i, j] = findmax(abs, C[i, j, :])[2]
        end
    end

    in_C = zeros(Float32, length(lon), length(lat))

    for i = 1:size(e_point_list)[1]
        for j = 1:size(i_point_list)[1]

            e_point = e_point_list[i, :]

            H(θ[j, i])
            in_C[e_point[1], e_point[2]] +=
                (C[j, i, θ[j, i]] * H(lags[θ[j, i]]) * (1 - (i == j)))
        end
    end

    fig = Figure()

    ga = GeoAxis(fig[1, 1]; coastlines = true)

    function pos(x)
        if x > 180
            return x - 360
        else
            return x
        end
    end

    plt = surface!(ga, lon, lat, in_C; shading = false, colormap = :inferno)

    cb = Colorbar(fig[1, 2], plt)

    close(generated_data)
    return C, lat, lon, i_point_list, e_point_list, in_C, fig, ga
end

C, lat, lon, i_point_list, e_point_list, in_C, fig, ga = main()
save("out.pdf", fig)
