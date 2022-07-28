#!/usr/bin/env julia

using Pkg
Pkg.add.(["HDF5","StatsBase","LoopVectorization"])

include("./nino.jl")
using .OurNino

function main()
    years=1977:1987
    data, lat, lon = get_data(; years=years)
    anomaly_data = get_anomaly(data; years=years)
	display(anomaly_data[1:100])

    enb_is, enb_js = find_inside_indeces(lat, lon, -170, -120, -5, 5)
	println()
    println(enb_is)
    println(enb_js)

    window = get_period(anomaly_data, 6, 1) 
    C, i_point_list, e_point_list = c_i_j(window, enb_is, enb_js)
	println(C[1:100])
	
	i_point_list = hcat(first.(i_point_list),last.(i_point_list))
	e_point_list = hcat(first.(e_point_list),last.(e_point_list))

	enb_is = Array(enb_is)
	enb_js = Array(enb_js)

	years = Array(years)

    in_C = in_weights(C, i_point_list, e_point_list, length(lon), length(lat))

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
    end

    return 0
end

main()
