using  DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics, LsqFit, LinearAlgebra

default(fontfamily = "Computer Modern",
background_color = :transparent,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c/spectrum_c"
cd(path)
sample, Nt= 100000000, 60.
fname = @sprintf "data_spectrum_sample%.1eNt%i.txt" sample Nt 

lines = readlines(fname)

p = plot(legend=:best, dpi=600)

t_dist = [parse(Int64, split(line, ' ')[1]) for line in lines[2:end-1]]
gap= [parse(Float64, split(line, ' ')[2]) for line in lines[2:end-1]]
gapv= [parse(Float64, split(line, ' ')[3]) for line in lines[2:end-1]]

scatter!(p, t_dist ,gap, yerror=gapv, markershape=:plus, label = L"N_t=N_s=60")
xlabel!(p,"spacing")
ylabel!(p,L"\frac{gap}{\hat{m}}")
title!("Plot of normalized Energy Gap vs Correlators Spacing")
display(p)
savefig(p, "gaps.png")

# Fit
p1 = plot(legend=:best, dpi=600)
model(x, p) = p[1] .+ 0 .* x
p0 = [0.0]
fit = curve_fit(model, t_dist, gap, gapv.^(-2), p0)
best_fit_params =coef(fit)
J = fit.jacobian
cov_matrix = inv(J' * J)
param_errors = sqrt.(Diagonal(cov_matrix))
# Stampa i risultati
println("Best fit Parameters ",best_fit_params )
println("Parameter errors: ", param_errors)


# Grafico dei dati e del fit
scatter!(p1, t_dist[1:10], gap[1:10], yerror= gapv[1:10],markershape=:plus,label ="data" )
plot!(p1, t_dist, model(t_dist, coef(fit)), label=L"y =  0.9951 $\pm$ 0.0004)", color=:violet)
xlabel!(p1,"spacing")
ylabel!(p1,L"\frac{gap}{\hat{m}}")
title!("Fit of normalized Energy Gap vs Correlators Spacing")
display(p1)
savefig(p1, "gaps_fit.png")