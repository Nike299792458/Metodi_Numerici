using CSV, DataFrames, DelimitedFiles, LaTeXStrings, Plots, Printf, Statistics, LsqFit, LinearAlgebra

default(
    fontfamily = "Computer Modern",
    background_color = :transparent,
    foreground_color = :black,
    background_color_legend = nothing,
    margin = 5Plots.mm
)

sample = 5000000

# T ≠ m
path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c"
cd(path)

plots = []
p2 = plot()
marker_shapes = [:cross, :cross, :cross]
time_division = [4,8,10]
ratio = 5

for (i, Nt) in enumerate(time_division)
    Nt_b = Nt * ratio
    fname = @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e.txt", Nt, Nt_b, sample)
    if isfile(fname)
        lines = readlines(fname)
        if length(lines) > 2
            Tonm = [parse(Float64, split(line, ',')[1]) for line in lines[2:end]]
            ϵ_norm = [parse(Float64, split(line, ',')[2]) for line in lines[2:end]]
            ϵ_normv = [parse(Float64, split(line, ',')[3]) for line in lines[2:end]]
            obs1_r = [parse(Float64, split(line, ',')[4]) for line in lines[2:end]]
            obs1v_r = [parse(Float64, split(line, ',')[5]) for line in lines[2:end]]

            p = plot()
            local y1 = 0.1866571
            scatter!(p2, Tonm, ϵ_norm, yerr = ϵ_normv, markershape = marker_shapes[i], label = "Nt=$Nt")
            scatter!(p, Tonm, obs1_r, yerr = obs1v_r, markershape = marker_shapes[i], label = "Nt=$Nt")
            plot!(p, [0.9, 1.1], [y1, y1], label = "theoretical value for T=m", lw = 2)
            push!(plots, p)
        else
            println("File $fname does not contain enough data.")
        end
    else
        println("File $fname not found.")
    end

end
# Combine all plots into a single figure with subplots
combined_plot = plot(plots..., layout = (length(plots), 1), size = (800, 1200))

for p in plots
   
    xlabel!("T/m")
    ylabel!(L"\frac{ϵ-p}{T^2}")
    title!("Temperature behavior")
    
end


display(combined_plot)

# Save combined plot
savefig(combined_plot, "o1suTsquared.png")

xlabel!(p2, "T/m")
ylabel!(p2, L"\frac{ϵ}{T^2}")
title!(p2, "Temperature behavior ")
y2 = π / 6
plot!(p2, [0, 2.5], [y2, y2], label = "high temerature limit", lw = 2)

display(p2)
savefig(p2, "ϵsuTsquared.png")



#T=m  
path = "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c/Tequalsm/"
cd(path)


p3 = plot(legend=:best, dpi=600)
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
    best_fit_params = coef(fit)
    lnsp=[0,0.04]
    plot!(p3, lnsp, model(lnsp, coef(fit)), label = "Linear fit ratio=$r", lw=2)
    J = fit.jacobian
    cov_matrix = inv(J' * J)
    param_errors = sqrt.(Diagonal(cov_matrix))
    println("Best fit parameters: ", best_fit_params)
    println("Parameter errors: ", param_errors)
    
    
end
title!(p3, L"$\frac{T}{m}=1$")
y1= 0.406123
plot!(p3, [0, 0.045], [y1, y1], label = "continuum limit", lw = 2)
display(p3) 
savefig(p3, "T=m.png")


p4 = plot(legend=:best, dpi=600)
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
    y1= 4*0.406123
    plot!(p4, [0, 0.045], [y1, y1], label = "continuum limit", lw = 2)
    title!(p4, L"$\frac{T}{m}=1$ wrong discretization")
    lnsp=[0,0.04]
    plot!(p4, lnsp, model(lnsp, coef(fit)), label = "Linear fit", lw=2)
    best_fit_params = coef(fit)
    J = fit.jacobian
    cov_matrix = inv(J' * J)
    param_errors = sqrt.(Diagonal(cov_matrix))
    println("Best fit parameters: ", best_fit_params)
    println("Parameter errors: ", param_errors)
    display(p4)
    savefig(p4, "doublers.png")
end


