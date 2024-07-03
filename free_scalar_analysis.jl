using ArgParse, DelimitedFiles, LinearAlgebra, Printf, Statistics
include("free_scalar.jl")


function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "sample"
            help = "number of steps in the simulation"
            required = true
            arg_type = Int
        "ratio"
            help = "ratio between spatial and temporal dimension"
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
        "doublers"
            help= "if we're simulating with wrong discretization"  
            default = false
            required = false
            arg_type = Int
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_c"])
            required = false
            arg_type = String
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    path = parsed_args["path"]
    doublers= parsed_args["doublers"]
    blocksize = parsed_args["blocksize"]
    therm = parsed_args["therm"]
    ratio = parsed_args["ratio"]
    sample = parsed_args["sample"]
    dfname = @sprintf("data_ratio=%i_sample=%.1e_doublers=%i.txt", ratio, sample, doublers)
    startp = @sprintf "free_scalar_th_sample=%.1eratio=%.i" sample ratio
    paths = filter(startswith(startp), readdir(path))
    temporal_dim = []
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
        ϵ_norm_j =(O1_j+O2_j-O3_j)/2

        O1_b_j = JackKnife(w[1][therm:end,4], blocksize)
        O2_b_j = JackKnife(w[1][therm:end,5], blocksize)
        O3_b_j = JackKnife(w[1][therm:end,6], blocksize)
        ϵ_norm_b_j =(O1_b_j+O2_b_j-O3_b_j)/2
 
        
        push!(temporal_dim, Nt)
        push!(O1, mean(O1_j))
        push!(ϵ_norm, (mean(ϵ_norm_j)- mean(ϵ_norm_b_j)))

        push!(O1v, std(O1_j, corrected = false).*sqrt(length(O1_j)-1))
        push!(ϵ_normv, std(ϵ_norm_j, corrected = false).*sqrt(length(ϵ_norm_j)-1))
    
    end

    w = open(joinpath([path, dfname]), "w") do io
        writedlm(io, ["temporal_dim" "O1" "O1v" "ϵ_norm" "ϵ_normv"], ",")
        writedlm(io, [temporal_dim O1 O1v ϵ_norm ϵ_normv], ",")
    end
    println("Done! Data stored in $(joinpath([path, dfname]))")
end
main()
