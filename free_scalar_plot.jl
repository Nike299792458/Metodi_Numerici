using Plots, LinearAlgebra, Statistics, Printf


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)
ratio= 4
sample= 100000
Nt=4
Nt_b=16
doublers= false
p1=plot()
p2=plot()



#T≠m
fname= @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
lines = readlines(fname)
Tonm = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
ϵ_norm= [parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
ϵ_normv= [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
x=Tonm
scatter!(p1, x, ϵ_norm, yerr=ϵ_normv, label="Dati ")
display(p1)




#=
#T=m  
temporal_division=[4,5,6,7,8,10]
ϵ_norm = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ϵ_norm_v = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x =[]
push!(x, temporal_division .^(-2))

for i= 1:length(temporal_division)
    println(i)
end
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
scatter!(p, x, ϵ_norm, yerr=ϵ_normv, label="Dati ")
display(p2)
=#
