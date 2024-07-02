using ArgParse, Dates, DelimitedFiles, Random, Printf, LinearAlgebra
include("free_scalar.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "ratio"
            help = "Ratio between spatial and temporal size of the lattice"
            required = true
            arg_type = Int
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        #=
        "mhat"
            help = "Mass in lattice units, a*m"
            required = true
            arg_type = Int
        =#
        "--path", "-p"
            help = "The path where files are stored"
            default = joinpath(["..", "simulations_c"])
            required = false
            arg_type = String
        
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    ratio = parsed_args["ratio"]
    sample = parsed_args["sample"]
    #mhat = parsed_args["mhat"]
    path = parsed_args["path"]
   

    
    println(@sprintf "Starting simulation: sample=%.1e ratio=%.i " sample ratio )

    timestamps=[4,5,6,7,8,10]
    calcolato=false
    for Nt in timestamps
        Ns=ratio*Nt
        # initializing...
        lattice = zeros(Float64, Nt,Ns)
        acc = 0
        
        
        # simulation parameters
        orsteps = 5
        measevery = 1
        Nt_bar= 10000 
        mhat=1/Nt
        stvol=Nt #stvol=Nt*Ns^{STDIM-1}
        stvol_bar= Ns*Nt_bar
        for i in 1:(STDIM-1)
            stvol=stvol*Ns
        end
        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf "free_scalar_th_sample=%.1eratio=%.iNt=%2.2iNs=%2.2i.txt" sample ratio Nt Ns
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1" "obs2" "obs3"], " ")
        end
        
        start = now()
        datafile = open(fr, "a")
        for iter in 1:sample
            for r in LinearIndices(lattice)
                acc+=heathbath!(lattice, r, mhat, Nt, Ns)
                for _ in 1:orsteps
                    acc+=overrelax!(lattice, r, mhat, Nt, Ns)
                end
            end

            if iter%measevery == 0
                obs1 = O1(stvol, mhat,lattice)-O1(stvol_bar, mhat,lattice)
                obs2 = O2(stvol, Nt,lattice)-O2(stvol_bar, Nt_bar,lattice)
                obs3 = O3(stvol,lattice)-O3(stvol_bar,lattice)
                writedlm(datafile, [obs1 obs2 obs3], " ")
            end
            
        end
        close(datafile)
        elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
        println("\n$(round(now(), Dates.Second));\nNₜ = $Nt,Ns = $Ns, elapsed time $(elapsed)\n")
    end
end
main()
