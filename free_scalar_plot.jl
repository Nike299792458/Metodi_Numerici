using CSV, DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)



sample= 5000000




doublers= false
#=

p1=plot()

#T≠m
# Definire un array di marker shapes
marker_shapes = [ :cross, :cross, :cross]

# Creiamo i plot principali
p1 = plot()
p2 = plot()

for (i, Nt) in enumerate(time_division)
    Nt_b = Nt * ratio
    fname = @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
    lines = readlines(fname)
    Tonm = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
    ϵ_norm = [parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
    ϵ_normv = [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
    obs1_r = [parse(Float64, split(line, ',')[4]) for line in lines[2:end]]
    obs1v_r = [parse(Float64, split(line, ',')[5]) for line in lines[2:end]]
    
    scatter!(p1, Tonm, ϵ_norm, yerr = ϵ_normv, markershape = marker_shapes[i], label = "Nt=$Nt")
    scatter!(p2, Tonm, obs1_r, yerr = obs1v_r, markershape = marker_shapes[i], label = "Nt=$Nt")
end

xlabel!(p1, "T/m")
ylabel!(p1, L"\frac{ϵ}{T^2}")
title!(p1, "Temperature behavior")
y1 = π / 6
plot!(p1, [0, 2.5], [y1, y1], label = "continuum limit", lw = 2)

xlabel!(p2, "T/m")
ylabel!(p2, L"\frac{ϵ-P}{T^2}")
title!(p2, "Temperature behavior")
y2 = 0.1866571
plot!(p2, [0, 2.5], [y2, y2], label = "continuum limit", lw = 2)

# Visualizziamo i plot finali
display(p1)
display(p2)
savefig(p1, "esuTsquared.png")
savefig(p2, "o1suTsquared.png")

=#

#T=m  
path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c/Tequalsm/"
cd(path)
p3=plot()


ratio=[4,5,6]
for (i,r) in enumerate(ratio)
    
    fname= @sprintf("Tequals_ratio=%.i_sample=%.1e.txt", r, sample)
    lines = readlines(fname)
    temporal_division=[5,6,8,10]
    #Nt=[parse(Int64, split(line, ',')[1]) for line in lines[2:end]]
    ϵ_norm=[parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
    ϵ_normv=[parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
    x=temporal_division .^(-2)
    scatter!(p3, x, ϵ_norm, yerr = ϵ_normv, markershape = :utriangle, label = "ratio=$r")  
end

display(p3)

p4=plot()


ratio=[5]
for (i,r) in enumerate(ratio)
    
    fname= @sprintf("Doublers_ratio=%.i_sample=%.1e.txt", r, sample)
    lines = readlines(fname)
    temporal_division=[5,6,8,10]
    #Nt=[parse(Int64, split(line, ',')[1]) for line in lines[2:end]]
    ϵ_norm=[parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
    ϵ_normv=[parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
    x=temporal_division .^(-2)
    scatter!(p4, x, ϵ_norm, yerr = ϵ_normv, markershape = :utriangle, label = "ratio=$r")  
end

display(p4)