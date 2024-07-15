using ArgParse, Dates, DelimitedFiles, Random, Printf, LinearAlgebra
include("free_scalar.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "ratio"
            help = "Ratio between spatial and temporal size of the lattice"
            required = true
            arg_type = Int
        "Nt"
            help = "temporal division"
            required = true
            arg_type = Int
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "--path", "-p"
            help = "The path where files are stored, if simulating T=m /../simulations_c/Tequalsm"
            default = joinpath(["..", "simulations_c"])
            required = false
            arg_type = String
        
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    ratio = parsed_args["ratio"]
    Nt = parsed_args["Nt"]
    sample = parsed_args["sample"]
    path = parsed_args["path"]

    #simulation parameters 
    Ns=ratio*Nt
    stvol=Nt #stvol=Nt*Ns^{STDIM-1}
    for i in 1:(STDIM-1)
         stvol=stvol*Ns
     end

    Nt_b = ratio*Nt
    Ns_b = Ns
    stvol_b = Nt_b # stvol=Nt*Ns^{STDIM-1}
    for i in 1:(STDIM-1)
        stvol_b = stvol_b * Ns_b
    end

    #generic
    orsteps = 5
    measevery = 5

    println(@sprintf "Starting simulation: sample=%.1e ratio=%.i Nt=%.i " sample ratio Nt )
    #Tonm =[1.0] #when simulating T=m
    #Tonm = collect(range(0.1, stop=2.5, length=16)) #Otherwise
    if Nt==8
        Tonm =[1.70]
    else
        Tonm=[1.54,2.50]
    end 

    for T_norm in Tonm
        # initializing...
        #=
        lattice = zeros(Float64, Nt ,Ns)
        acc=0
        mhat=1/(Nt*T_norm)
        
        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf "fs_th_sample=%.1eratio=%.iNt=%2.2iTonm=%2.2f.txt"  sample ratio Nt T_norm
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
                obs1 = O1(stvol, mhat,lattice)
                obs2 = O2(stvol, Nt, Ns, lattice)
                obs3 = O3(stvol, Nt, lattice)
                writedlm(datafile, [obs1 obs2 obs3 ], " ")
            end
            
        end
        close(datafile)
     =#
        lattice_b = zeros(Float64, Nt_b, Ns_b)
        acc = 0

        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf("fs_th_sample=%.1eratio=%.iNt_b=%2.2iNt=%2.2iTonm=%2.2f.txt" , sample, ratio, Nt_b, Nt, T_norm)
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1_b" "obs2_b" "obs3_b"], " ")
        end

        datafile = open(fr, "a")
        for iter in 1:sample
            for r in LinearIndices(lattice_b)
                acc += heathbath!(lattice_b, r, mhat, Nt_b, Ns_b)
                for _ in 1:orsteps
                    acc += overrelax!(lattice_b, r, mhat, Nt_b, Ns_b)
                end
            end

            if iter % measevery == 0
                obs1_b = O1(stvol_b, mhat, lattice_b)
                obs2_b = O2(stvol_b, Nt_b, Ns_b, lattice_b)
                obs3_b = O3(stvol_b, Nt_b, lattice_b)
                writedlm(datafile, [obs1_b obs2_b obs3_b ], " ")
            end
        end
        close(datafile)

        elapsed = Dates.canonicalize(Dates.round((now() - start), Dates.Second))
        println("\n$(round(now(), Dates.Second));\nNₜ = $Nt,Tonm = $T_norm, elapsed time $(elapsed)\n")
    end
    
end
main()


