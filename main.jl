#!/usr/bin/env julia

using Pkg
Pkg.add.(["HDF5","StatsBase"])

include("./nino.jl")
using .OurNino

function main()
    data, lat, lon = get_data()
    anomaly_data = get_anomaly(data)
	display(anomaly_data[1:100])

    enb_is, enb_js = find_inside_indeces(lat, lon, 120, 170, -5, 5)
	println()
    println(enb_is)
    println(enb_js)

    window = get_period(data, 35, 1) 
    println(c_i_j(window, enb_is, enb_js)[1,2])

    return 0
end

main()
