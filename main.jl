#!/usr/bin/env julia

include("./nino.jl")

function main()
    anomaly_data=get_anomaly()
    display(anomaly_data)
    return 0
end

main()
