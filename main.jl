#!/usr/bin/env julia

include("./nino.jl")

function main()
    anomaly_data=get_anomaly()
	display(anomaly_data[1:1000])
    return 0
end

main()
