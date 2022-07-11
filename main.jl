#!/usr/bin/env julia

include("./nino.jl")

function main()
    data, lat, lon = get_data()
    anomaly_data = get_anomaly(data)
	display(anomaly_data[1:1000])

    enb_is, enb_js = find_inside_indeces(lat, lon, 120, 170, 5, -5)
    display(enb_is, enb_js)

    

    return 0
end

main()
