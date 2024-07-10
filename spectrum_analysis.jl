using ArgParse, DelimitedFiles, LinearAlgebra, Printf, Statistics, Plots, DataFrames
include("spectrum_c.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "blocksize"
            help = "the number of points in a block"
            required = true
            arg_type = Int
        "therm"
            help = "the number of points to discard as thermalization"
            required = true
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
    blocksize = parsed_args["blocksize"]
    therm = parsed_args["therm"]
    sample = parsed_args["sample"]
    startp = @sprintf "spectrum_sample=%.1e" sample
    paths = filter(startswith(startp), readdir(path))
    Tonm=5
    Nt=60
    Ns=10
    mhat=1/(Nt*Tonm)

    for fname in paths
        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end

        spacing=Vector{Int}()
        for el in w[2][2:end]
            append!(spacing, parse(Int, el[4:end]))
        end

        # Estrazione della matrice cormat
        cormat = w[1][therm+1:end,1:end]#serve il -2 altrimenti prende uno spazio, che nel file però non c'è 
        
        # Creazione dell'array datajack con le dimensioni adeguate
        datajack = zeros(size(cormat, 1) ÷ blocksize, size(cormat, 2))
        println(typeof(datajack))

        # Popolamento di datajack usando la funzione JackKnife
        for i in 1:size(cormat, 2)
            datajack[:, i] = JackKnife(cormat[:, i], blocksize)
        end
       
        gaps = log.(datajack./circshift(datajack, (0,-4)))
        errs = std(gaps, dims = 1, corrected = false).*sqrt(size(datajack,1)-1)
        gaps = mean(gaps, dims = 1)
       
        data = DataFrame(
            sp = spacing[1:end],
            gap = gaps[1:end]./mhat,
            err = errs[1:end]./mhat,
        )
        if !isdir(joinpath([path, "data"]))
            mkpath(joinpath([path, "data"]))
        end
                
        f = @sprintf "data_spectrum_sample%.1eNt%i.txt" sample Nt 
        touch(joinpath([path, f]))        
        w = open(joinpath([path, f]), "w") do io
        writedlm(io, ["spacing"  "gaps" "errs" ], " ")
        writedlm(io, collect(eachrow(data)), " ")
        end
    end
end    

main()

