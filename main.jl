#!/usr/bin/env julia

using Pkg
Pkg.add.(["HDF5","StatsBase","LoopVectorization"])

include("./nino.jl")
using .OurNino

function main()
    years=1972:1992
    data, lat, lon = get_data(; years=years)
    anomaly_data = get_anomaly(data; years=years)
	display(anomaly_data[1:100])

    enb_is, enb_js = find_inside_indeces(lat, lon, 120, 170, -5, 5)
	println()
    println(enb_is)
    println(enb_js)

    window = get_period(data, 11, 1) 
    C, i_point_list, e_point_list = c_i_j(window, enb_is, enb_js)
	println(C[1:100])

    return 0
end

main()
