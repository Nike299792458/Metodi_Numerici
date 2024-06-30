using Plots, LinearAlgebra, Statistics


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)
sample = 100000
ratio=8
fname = @sprintf "data_ratio=4" 


lines = readlines("data_ratio=4")
lines=lines[2:end]
Nt = [parse(Float64, split(line, ',')[1]) for line in lines]
ϵ_norm= [parse(Float64, split(line, ',')[4]) for line in lines]
ϵ_normv= [parse(Float64, split(line, ',')[5]) for line in lines]


x=Nt .^(-2)
scatter(x, ϵ_norm, yerr =ϵ_normv, title="Scatter Plot", xlabel="x", ylabel="y", label="Dati")


