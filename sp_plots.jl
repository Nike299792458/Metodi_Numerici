using DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics


default(fontfamily = "Computer Modern",
background_color = :transparent,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)

path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
sample, Nt= 1000000, 20.
fname = @sprintf "data_spectrum_sample%.1eNt%i.txt" sample Nt 

data = readdlm(joinpath(path, fname))
df = DataFrame(data, :auto)

p = plot(legend=:best, dpi=600)

scatter!(p, df[!,:sp], df[!,:gaps], yerr = df[!,:errs], markershape=:plus, label = L"x-x")

xlabel!("Spacing")
ylabel!(L"\Delta E/(\hbar \omega)")
ylims!(0.9, 4)

title!("Plot of Energy Gaps vs Correlators Spacing")

display(p)
# savefig(p, "..\\imgs_b\\gaps.png")