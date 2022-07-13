#!/usr/bin/env julia

using HDF5, GLMakie, GeoMakie

function H(x)
    if x >= 0
        return x
    else
        return typeof(x)(0)
    end
end

function main()
    generated_data = h5open("generated_data.h5", "r") 

    C = read(generated_data["C"])
    size_C = size(C) 

    lat = read(generated_data["lat"])
    lon = read(generated_data["lon"])

    e_point_list = read(generated_data["e_point_list"])
    i_point_list = read(generated_data["i_point_list"])

    θ = Array{Float32}(undef, size_C[1], size_C[2])
    for i in 1:size_C[1]
        for j in 1:size_C[2]
            θ[i,j] = C[i,j,findmax(abs,C[i,j,:])[2]]
        end
    end
    
    #in_C = Array{Float32}(undef, length(lat), length(lon))
    in_C = zeros(Float32, length(lon), length(lat))

    for i in 1:size(e_point_list)[1]
        for j in 1:size(i_point_list)[1]

            #i_point = i_point_list[j,:]
            #in_C[i_point[1], i_point[2]]=100
            
            e_point = e_point_list[i,:]

            H(θ[j,i])
            in_C[e_point[1], e_point[2]]
            in_C[e_point[1], e_point[2]] += θ[j,i] 
        end
    end
    

    fig = Figure()

    ga = GeoAxis(
		fig[1,1];
		coastlines = true,
	)

    plt = surface!(ga, lon, lat, in_C; shading=false, colormap=:inferno)
    
    #anomaly_data = read(generated_data["anomaly_data"])
    #plt = surface!(ga, lon, lat, anomaly_data[:,:,1,11]; shading=false, colormap=:inferno)
    
    #plt2 = plot!(ga,Point2f(-170,-5),Point2f(-120,5))
    #xlims!(ga,-90,90)
    cb = Colorbar(fig[1,2], plt)

    close(generated_data)
    return C, lat, lon, i_point_list, e_point_list, in_C, fig, ga
end

C, lat, lon, i_point_list, e_point_list, in_C, fig, ga = main()
fig