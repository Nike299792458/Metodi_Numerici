using ArgParse, DelimitedFiles, LinearAlgebra, Plots, Printf, Statistics
include("free_scalar.jl")


function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "therm"
            help = "the number of points to discard as thermalization"
            required = true
            arg_type = Int
        "blocksize"
            help = "the number of points in a block"
            required = true
            arg_type = Int
        "Nt"
            help = "Temporal dimension of the lattice"
            required = true
            arg_type = Int       
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_c"])
            required = false
            arg_type = String
        "--name", "-n"
            help = "the name of data file"
            default = "data.txt"
            required = false
            arg_type = String
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    path = parsed_args["path"]
    blocksize = parsed_args["blocksize"]
    therm = parsed_args["therm"]
    sample = parsed_args["sample"]
    dfname = parsed_args["name"]
    startp = @sprintf "free_scalar_th_sample=%.1e" sample
    paths = filter(startswith(startp), readdir(path))
    size_list = []
    O1 = []
    ϵ_norm = []
    O1v = []
    ϵ_normv = []
    for (i,fname) in enumerate(paths)
        Nt = parse(Int, fname[end-10:end-9])

        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end

        O1_j = JackKnife(w[1][therm:end,1], blocksize)
        O2_j = JackKnife(w[1][therm:end,2], blocksize)
        O3_j = JackKnife(w[1][therm:end,3], blocksize)
        ϵ_norm_j = (O1_j+O2_j-O3_j)/2
        
        
        push!(size_list, Nt)
        push!(O1, mean(O1_j))
        push!(ϵ_norm, mean(ϵ_norm_j))

        push!(O1v, std(O1_j, corrected = false).*sqrt(length(O1_j))-1)
        push!(ϵ_normv, std(ϵ_norm_j, corrected = false).*sqrt(length(O1_j)-1))
    
    end

    w = open(joinpath([path, dfname]), "w") do io
        writedlm(io, ["size_list" "O1" "O1v" "ϵ_norm" "ϵ_normv"], ",")
        writedlm(io, [size_list O1 O1v ϵ_norm ϵ_normv], ",")
    end
    println("Done! Data stored in $(joinpath([path, dfname]))")
end

main()
#=
stampa qualcosa ma c'è da capire cosa 
e perché sono 4 righe non so come fa le medie=#