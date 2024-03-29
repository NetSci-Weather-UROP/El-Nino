#!/usr/bin/env julia

if !haskey(ENV, "PKG_UPD")
    using Pkg
    Pkg.add.(["HDF5", "StatsBase", "LoopVectorization"])
end

include("./nino.jl")
using .OurNino

if haskey(ENV, "ANOM_IMPL")
    const anom_impl = ENV["ANOM_IMPL"]
else
    const anom_impl = "1"
end

if haskey(ENV, "C_IMPL")
    const c_impl = ENV["C_IMPL"]
else
    const c_impl = "1"
end

function main()
    if anom_impl == "2"
        years = 1948:2021
    else
        years = 1967:1977
    end

    data, lat, lon = get_data(; years = years)

    if anom_impl == "1"
        anomaly_data = get_anomaly(data; years = years)
    elseif anom_impl == "2"
        anomaly_data = get_anomaly2(data; years = years)
    end

    display(anomaly_data[1:100])

    enb_is, enb_js = find_inside_indeces(lat, lon, -170, -120, -5, 5)
    println()
    println(enb_is)
    println(enb_js)

    window = get_period(anomaly_data, findfirst(isequal(1971), years), 348)

    if c_impl == "1"
        C, i_point_list, e_point_list = c_i_j(window, enb_is, enb_js)
    elseif c_impl == "2"
        C, i_point_list, e_point_list = c_i_j2(window, enb_is, enb_js)
    elseif c_impl == "3"
        C, i_point_list, e_point_list = c_i_j3(window, enb_is, enb_js)
    elseif c_impl == "4"
        C, i_point_list, e_point_list = c_i_j4(window, enb_is, enb_js)
    end

    println(C[1:100])

    i_point_list = hcat(first.(i_point_list), last.(i_point_list))
    e_point_list = hcat(first.(e_point_list), last.(e_point_list))

    enb_is = Array(enb_is)
    enb_js = Array(enb_js)

    years = Array(years)

    θ = get_θ(C, size(C))
    in_C = in_weights(C, i_point_list, e_point_list, length(lon), length(lat), θ)

    OurNino.HDF5.h5open("generated_data.h5", "w") do file
        OurNino.HDF5.@write file C
        OurNino.HDF5.@write file i_point_list
        OurNino.HDF5.@write file e_point_list
        OurNino.HDF5.@write file enb_is
        OurNino.HDF5.@write file enb_js
        OurNino.HDF5.@write file anomaly_data
        OurNino.HDF5.@write file lat
        OurNino.HDF5.@write file lon
        OurNino.HDF5.@write file years
        OurNino.HDF5.@write file in_C
        OurNino.HDF5.@write file window
    end

    return 0
end

main()
