using CSV, DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)
ratio= 5
sample= 5000000
Nt=4
Nt_b=20
doublers= false
p1=plot()
p2=plot()
p3=plot()
p4=plot()



#T≠m
fname= @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
lines = readlines(fname)
Tonm = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
ϵ_norm= [parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
ϵ_normv= [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
x=Tonm
scatter!(p1,ϵ_norm ,Tonm, yerr = ϵ_normv, markershape=:plus, label = "Nt=4")
xlabel!(p1,"T/m")
ylabel!(p1,L"\frac{ϵ}{T^2}")
#ylims!(0.9, 4)
#annotate!(p, 2.8, 3.33, text(L"Where $A = x^3-\frac{3}{2}x$", :black,:left, 8))
title!("Termodynamics")
display(p1)

fname= @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
lines = readlines(fname)
Tonm = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
obs1_r= [parse(Float64, split(line, ',')[4]) for line in lines[2:end]]
obs1v_r= [parse(Float64, split(line, ',')[5]) for line in lines[2:end]]
x=Tonm
xlabel!(p2,"T/m")
ylabel!(p2,L"\frac{ϵ-P}{T^2}")
scatter!(p2,x ,obs1_r, yerr = obs1v_r, markershape=:plus, label = "Nt=4")
title!("Termodynamics")
display(p2)


#T=m  
#=
temporal_division=[5,6,8,10]
ϵ_norm = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ϵ_norm_v = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x =[]
push!(x, temporal_division .^(-2))

for i= 1:length(temporal_division)
    local Nt = temporal_division[i]
    local Nt_b= Nt*ratio
    local fname= @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
    local lines = readlines(fname)
    line = lines[2]
    ϵ_norm[i]= parse(Float64, split(line, ',')[2])  
    ϵ_normv[i]= parse(Float64, split(line, ',')[3]) 
end


scatter!(p3, x, ϵ_norm, yerr=ϵ_normv, label="Dati ")
display(p3)
#doublers
for i= 1:length(temporal_division)
    local Nt = temporal_division[i]
    local Nt_b= Nt*ratio
    local fname= @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
    local lines = readlines(fname)
    line = lines[2]
    ϵ_norm[i]= parse(Float64, split(line, ',')[2])  
    ϵ_normv[i]= parse(Float64, split(line, ',')[3]) 
end

println(ϵ_norm)
scatter!(p4, x, ϵ_norm, yerr=ϵ_normv, label="Dati ")
display(p4)

=#