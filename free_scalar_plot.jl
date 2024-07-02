using Plots, LinearAlgebra, Statistics, Printf


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)
ratio= [4]
p=plot()
for r in ratio
    println(r)
    fname= @sprintf("data_ratio=%i", r)
    lines = readlines(fname)
    Nt = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
    ϵ_norm= [parse(Float64, split(line, ',')[4]) for line in lines[2:end]]
    ϵ_normv= [parse(Float64, split(line, ',')[5]) for line in lines[2:end]]
    x=Nt .^(-2)
    scatter!(p, x, 2 .* ϵ_norm, yerr=ϵ_normv, label="Dati $r")
    display(p)
end


