using ArgParse, Dates, DelimitedFiles, Random, Printf, LinearAlgebra
include("doublers.jl")

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
    Nt = parsed_args["Nt"]
    sample = parsed_args["sample"]
    path = parsed_args["path"]

    #simulation parameters (T/m independant)
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
    measevery = 10

    println(@sprintf "Starting simulation: sample=%.1e ratio=%.i Nt=%.i " sample ratio Nt )
    Tonm =[1] #when simulating T=m
    for T_norm in Tonm
        # initializing...
        lattice = zeros(Float64, Nt ,Ns)
        acc=0
        mhat=1/(Nt*T_norm)
        
        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf "doublers_sample=%.1eratio=%.iNt=%2.2iTonm=%2.2f.txt"  sample ratio Nt T_norm
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1d" "obs2d" "obs3d"], " ")
        end
        
        start = now()
        datafile = open(fr, "a")
        for iter in 1:sample
            for r in LinearIndices(lattice)
                acc+=heathbathd!(lattice, r, mhat, Nt, Ns)
                for _ in 1:orsteps
                    acc+=overrelaxd!(lattice, r, mhat, Nt, Ns)
                end
            end

            if iter%measevery == 0
                obs1d = O1d(stvol, mhat,lattice)
                obs2d = O2d(stvol, Nt, Ns, lattice)
                obs3d = O3d(stvol, Nt, lattice)
                writedlm(datafile, [obs1d obs2d obs3d ], " ")
            end
            
        end
        close(datafile)

        lattice_b = zeros(Float64, Nt_b, Ns_b)
        acc = 0

        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf("doublers_sample=%.1eratio=%.iNt_b=%2.2iNt=%2.2iTonm=%2.2f.txt" , sample, ratio, Nt_b, Nt, T_norm)
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1d_b" "obs2d_b" "obs3d_b"], " ")
        end

        datafile = open(fr, "a")
        for iter in 1:sample
            for r in LinearIndices(lattice_b)
                acc += heathbathd!(lattice_b, r, mhat, Nt_b, Ns_b)
                for _ in 1:orsteps
                    acc += overrelaxd!(lattice_b, r, mhat, Nt_b, Ns_b)
                end
            end

            if iter % measevery == 0
                obs1d_b = O1d(stvol_b, mhat, lattice_b)
                obs2d_b = O2d(stvol_b, Nt_b, Ns_b, lattice_b)
                obs3d_b = O3d(stvol_b, Nt_b, lattice_b)
                writedlm(datafile, [obs1d_b obs2d_b obs3d_b], " ")
            end
        end
        close(datafile)

        elapsed = Dates.canonicalize(Dates.round((now() - start), Dates.Second))
        println("\n$(round(now(), Dates.Second));\nNₜ = $Nt,Tonm = $T_norm, elapsed time $(elapsed)\n")
    end
    
end
main()


