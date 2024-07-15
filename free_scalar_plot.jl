using CSV, DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics, LsqFit, LinearAlgebra


default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)
sample=5000000

#T≠m 
path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)

p1 = plot()
p2 = plot()
marker_shapes = [ :cross, :dot, :utriangle]
time_division=[4,8,10]
ratio=5

for (i, Nt) in enumerate(time_division)
    Nt_b = Nt * ratio
    fname = @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e.txt", Nt, Nt_b, sample)
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
ylabel!(p2, L"\frac{ϵ-p}{T^2}")
title!(p2, "Temperature behavior")
y2 = 0.1866571
plot!(p2, [0, 2.5], [y2, y2], label = "continuum limit", lw = 2)


display(p1)
display(p2)
savefig(p1, "esuTsquared.png")
savefig(p2, "o1suTsquared.png")

#=
#T=m  
path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c/Tequalsm/"
cd(path)


p3=plot()
ratio=[4,5,6]
for (i,r) in enumerate(ratio)
    fname= @sprintf("Tequals_ratio=%.i_sample=%.1e.txt", r, sample)
    temporal_division=[5,6,8,10]
    lines = readlines(fname)
    ϵ_norm=[parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
    ϵ_normv=[parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
    
    x=temporal_division .^(-2)
    

    #fit
    model(x, p) = p[1] .+ p[2] .* x
    p0 = [0.0, 0.0]  
    fit = curve_fit(model, x, ϵ_norm, ϵ_normv.^(-2), p0)  
    scatter!(p3, x, ϵ_norm, yerr = ϵ_normv, markershape = :utriangle, label = "ratio=$r", xlabel = L"$\frac{1}{Nt^2}$", ylabel = L"$ϵ/T^2$")
    title!(p3, L"$\frac{T}{m}=1$")
    plot!(p3, x, model(x, coef(fit)), label = "Linear fit ratio=$r", lw=2)
    best_fit_params = coef(fit)
    J = fit.jacobian
    cov_matrix = inv(J' * J)
    param_errors = sqrt.(Diagonal(cov_matrix))
    println("Best fit parameters: ", best_fit_params)
    println("Parameter errors: ", param_errors)
    
    
end
display(p3) 
savefig(p3, "T=m.png")


p4=plot()
ratio=[5]
for (i,r) in enumerate(ratio)
    
    fname= @sprintf("data_doublers_sample=%.1eratio=%.i.txt", sample, r)
    lines = readlines(fname)
    temporal_division=[6,8,10]
    #Nt=[parse(Int64, split(line, ',')[1]) for line in lines[2:end]]
    ϵ_norm=[parse(Float64, split(line, ',')[2]) for line in lines[3:end]]
    ϵ_normv=[parse(Float64, split(line, ',')[3]) for line in lines[3:end]]
    x=temporal_division .^(-2)
     
    
    #fit
    model(x, p) = p[1] .+ p[2] .* x
    p0 = [0.0, 0.0]  
    fit = curve_fit(model, x, ϵ_norm, ϵ_normv.^(-2), p0)  
    intercept, slope = coef(fit)
  
    p4 = scatter(x, ϵ_norm, yerr = ϵ_normv, markershape = :plus, label = "ratio=$r", xlabel = L"$\frac{1}{Nt^2}$", ylabel = L"$ϵ/T^2$")
    title!(p4, L"$\frac{T}{m}=1$ wrong discretization")
    plot!(p4, x, model(x, coef(fit)), label = "Linear fit", lw=2)
    best_fit_params = coef(fit)
    J = fit.jacobian
    cov_matrix = inv(J' * J)
    param_errors = sqrt.(Diagonal(cov_matrix))
    println("Best fit parameters: ", best_fit_params)
    println("Parameter errors: ", param_errors)
    display(p4)
    savefig(p4, "doublers.png")
end

=#
