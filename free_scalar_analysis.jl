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
        "Nt"
            help = "tempomporal division"
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
            arg_type = Bool
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath([ "..", "simulations_c"])
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
    Nt = parsed_args["Nt"]
    Nt_b=Nt*ratio
    
    dfname = @sprintf("data_Nt=%2.2i_sample=%.1e_doublers=%i.txt", Nt, sample, doublers)
    startp = @sprintf "fs_th_sample=%.1eratio=%.iNt=%2.2iTonm=" sample ratio Nt  
    paths = filter(startswith(startp), readdir(path))
    Tonm = []
    sum_obs_j=[]
    sum_obs_b_j=[]
    ϵ_norm = []
    ϵ_norm_b =[]
    ϵ_norm_r = []
    ϵ_normv = []
    ϵ_normv_b =[]
    ϵ_normv_r = []
    obs1 = []
    obs1v = []
    obs1_b = []
    obs1v_b = []
    obs1_r = []
    obs1v_r = []




    for (i,fname) in enumerate(paths)
        T_norm =  parse(Float64, fname[end-7:end-4])

        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end

        O1_j = JackKnife(w[1][therm:end,1], blocksize)
        O2_j = JackKnife(w[1][therm:end,2], blocksize)
        O3_j = JackKnife(w[1][therm:end,3], blocksize)
        sum_obs_j =(O1_j+O2_j-O3_j)
        
        

        push!(Tonm , T_norm)
        push!(ϵ_norm , Nt*Nt*mean(sum_obs_j)/2)
        push!(ϵ_normv , Nt*Nt*std(sum_obs_j, corrected = false).*sqrt(length(sum_obs_j)-1))
        push!(obs1 , Nt*Nt*mean(O1_j))
        push!(obs1v , Nt*Nt*std(O1_j, corrected = false).*sqrt(length(O1_j)-1))

    end

    dfname = @sprintf("data_Nt=%2.2i_Nt_b=%2.2i_sample=%.1e_doublers=%i.txt", Nt, Nt_b, sample, doublers)
    startp = @sprintf("fs_th_sample=%.1eratio=%.iNt_b=%2.2iNt=%2.2iTonm=" ,sample, ratio, Nt_b, Nt) 
    paths = filter(startswith(startp), readdir(path))
    for (i,fname) in enumerate(paths)
        T_norm= parse(Float64, fname[end-7:end-4]) 

        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end

        O1_b_j = JackKnife(w[1][therm:end,1], blocksize)
        O2_b_j = JackKnife(w[1][therm:end,2], blocksize)
        O3_b_j = JackKnife(w[1][therm:end,3], blocksize)
        sum_obs_b_j =(O1_b_j+O2_b_j-O3_b_j)
        

        
        push!(ϵ_norm_b, Nt*Nt*mean(sum_obs_b_j)/2)
        push!(ϵ_normv_b, Nt*Nt*std(sum_obs_b_j, corrected = false).*sqrt(length(sum_obs_b_j)-1))
        push!(obs1_b , Nt*Nt*mean(O1_b_j))
        push!(obs1v_b , Nt*Nt*std(O1_b_j, corrected = false).*sqrt(length(O1_b_j)-1))

    end
    ϵ_norm_r= ϵ_norm - ϵ_norm_b
    ϵ_normv_r = sqrt.((ϵ_normv ./ ϵ_norm).^2 .+ (ϵ_normv_b ./ ϵ_norm_b).^2) .* ϵ_norm_r
    obs1_r= obs1 - obs1_b
    obs1v_r = sqrt.((obs1v ./ obs1).^2 .+ (obs1v_b ./ obs1_b).^2) .* obs1_r
    println(ϵ_norm_r)
    println(ϵ_normv_r)


    w = open(joinpath([path, dfname]), "w") do io
        writedlm(io, ["Tonm"  "ϵ_norm" "ϵ_normv" "obs1" "obs1v" ], ",")
        writedlm(io, [Tonm ϵ_norm_r ϵ_normv_r obs1_r obs1v_r] , ",")
    end
    println("Done! Data stored in $(joinpath([path, dfname]))")
end
main()


#=
Nt^(STDIM)*O1
=#