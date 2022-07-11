#!/usr/bin/env julia

include("./nino.jl")

function main()
    data = get_data()
    anomaly_data, lat, lon=get_anomaly(data)
	display(anomaly_data[1:1000])

    enb_i, enb_j = find_inside_indeces(lat, lon, 120, 170, 5, -5)
    display(enb_i, enb_j)
    
    return 0
end

main()
