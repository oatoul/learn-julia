using Base: String
using StatsBase: push!, length
using CartesianGeneticProgramming
using Cambrian
using CSV
using DataFrames
using StatsBase
import Cambrian.mutate


function evaluate!(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    accuracy = 0.0
    for i in 1:size(X, 1)
        out = process(ind, collect(X[i]))
        if out[1] == Y[i]
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end


function max_selection(pop::Array{<:Individual})
    sort(unique(pop), by=i ->(i.fitness, -i.sparsity, -i.n_active))[end]
end


function getBN!(df::DataFrame, fitness::Float64, exptect_conn::Set, universe_conn::Set)
    BN_mlp = BooleanNetwork(exptect_conn=exptect_conn, universe_conn=universe_conn)
    BN_std = BooleanNetwork(exptect_conn=exptect_conn, universe_conn=universe_conn)
    BN_fit = BooleanNetwork(exptect_conn=exptect_conn, universe_conn=universe_conn)

    ndf = names(df)
    l_idx = 3
    h_idx = size(ndf)[1]
    
    X = copy(Tuple.(eachrow(df)))
    pop!(X)

    for i = l_idx:h_idx
        target = ndf[i]
        println("Calculate $(i): $(target)")

        "data setup for target gene"
        Y = copy(df[!,target])
        popfirst!(Y)

        
        fit(ind::CGPInd) = evaluate!(ind, X, Y)
        e = CGPEvolution_archive(cfg, fit)

        step_archive!(e)
        while(e.population[end].fitness[1] < fitness && e.gen < e.config.n_gen)
            step_archive!(e)
        end

        println("Fitness: $(e.elites[end].fitness[1])")

        # normal selection
        ind_std = e.population[end]
        
        # archive fit selection
        if size(e.elites)[1] > 1
            ind_fit = max_selection(e.elites)
        else
            ind_fit = e.elites[end]
        end
        
        # archive mlp selection
        if size(e.elites)[1] > 1
            ind_mlp = get_best_from_archive!(e.elites, l_idx, h_idx, i)
        else
            ind_mlp = e.elites[end]
        end
        
        set_std = get_active_connections!(ind_std, l_idx, h_idx)
        set_fit = get_active_connections!(ind_fit, l_idx, h_idx)
        set_mlp = get_active_connections!(ind_mlp, l_idx, h_idx)

        add_dynamic_consistency!(BN_std, ind_std.fitness[1])
        add_dynamic_consistency!(BN_fit, ind_fit.fitness[1])
        add_dynamic_consistency!(BN_mlp, ind_mlp.fitness[1])

        store_BN(set_std, BN_std, ndf, target)
        store_BN(set_fit, BN_fit, ndf, target)
        store_BN(set_mlp, BN_mlp, ndf, target)
    end

    println(BN_std.dynamic_consistency)
    println(BN_fit.dynamic_consistency)
    println(BN_mlp.dynamic_consistency)

    BN_std, BN_fit, BN_mlp
end


function store_BN(set::Set, BN::BooleanNetwork, ndf::Array, target::String)
    for k in set
        conn = ndf[k] * "_" * target
        push!(BN.actual_conn, conn)
    end
    BN
end


function get_connections_by_vote!(vote_ratio::Float64, elites::Array{CGPInd}, low::Int64, high::Int64)
    res_set = Set()
    res_all = []
    min_vote = length(elites) * vote_ratio
    println("min vote is $(min_vote)")
    for e in elites
        r = collect(get_active_connections!(e, low, high))
        res_all = vcat(res_all, r)
    end
    count = StatsBase.countmap(res_all)
    for (k,v) in count
        if v >= min_vote
            push!(res_set,k)
        end
    end
    res_set
end

function get_best_from_archive!(elites::Array{CGPInd}, low::Int64, high::Int64, target_idx::Int64)
    res_all = []
    shift = low - 1
    for e in elites
        loss = evaluate_mlp(e, low, high, target_idx-shift)
        res_all = vcat(res_all, (loss, e))
        # println("Fitness $(e.fitness) sparsity $(e.sparsity) active $(e.n_active) loss: $(loss)")
    end

    best = sort(res_all, by=((x,y),)->(-x,))[end]
    loss = best[1]
    best = best[2]
    println("Best: Fitness $(best.fitness) sparsity $(best.sparsity) active $(best.n_active) loss $(loss)")

    best
end

function evaluate_mlp(ind::CGPInd, low::Int64, high::Int64, target_idx::Int64)
    in_idx = collect(Int64, get_active_connections_no_shift!(ind, low, high))
    train_mlp_flux(df_float, in_idx, target_idx)
end


"
(TP + TN) / (TP + FP + FN + TN)

TP (true positive) correctly predicted connections
FP (false positive) incorrectly predicted connections
FN (false negative) non-inferred connections
TN (true negative) correct negative predictions

"
function get_structural_accuracy(expect::Set, actual::Set, universe::Set)
    TP = length(intersect(actual, expect))
    FP = length(setdiff(actual, expect))
    FN = length(setdiff(expect, actual))
    expect_negative = setdiff(universe, expect)
    actual_negative = setdiff(universe, actual)
    TN = length(intersect(actual_negative, expect_negative))

    acc = (TP + TN) / (TP + FP + FN + TN)
    println("TP : $(TP)")
    println("FP : $(FP)")
    println("FN : $(FN)")
    println("TN : $(TN)")
    println("Structural accuracy : $(acc)")

    acc
end

function get_universe_set(inputs::Any)
    res = Set()
    for i in inputs
        for j in inputs
            re = i * "_" * j
            push!(res, re)
        end
    end
    res
end

function get_expect_Silico2()
    res = Set()
    push!(res, "G1_G2")

    push!(res, "G2_G4")
    push!(res, "G2_G5")
    push!(res, "G2_G6")
    push!(res, "G2_G7")
    push!(res, "G2_G8")

    push!(res, "G3_G4")
    push!(res, "G3_G5")

    push!(res, "G4_G3")
    push!(res, "G4_G5")

    push!(res, "G5_G2")
    push!(res, "G5_G3")

    push!(res, "G6_G7")

    push!(res, "G7_G6")

    push!(res, "G9_G5")

    push!(res, "G10_G3")

    res
end


cfg = get_config("cfg/Silico2.yaml")
println(cfg)
mutate(ind::CGPInd) = goldman_mutate(cfg, ind)

df_origin_bool = DataFrame(CSV.File("data/Silico2_bool.tsv",drop=["Time"],type=Bool))
df_origin_int = DataFrame(CSV.File("data/Silico2_bool.tsv",drop=["Time"],type=Int))

df_bool = copy(df_origin_bool)
insertcols!(df_bool, 1, :T0 => false, :T1 => true)

println(df_bool)

df_origin_float = DataFrame(CSV.File("data/Silico2.tsv",drop=["Time"],type=Float32))
df_float = Matrix(df_origin_float)'

universe = get_universe_set(names(df_origin_bool))

expect = get_expect_Silico2()

# for i = 1:10
#     println("########## Starting iteration $(i) #################")

#     BN_std, BN_fit, BN_mlp = getBN!(df_bool, 0.99, expect, universe)

#     println("Structural acc for STD")
#     evaluate_bn!(BN_std)

#     println("Structural acc for FIT")
#     evaluate_bn!(BN_fit)

#     println("Structural acc for MLP")
#     evaluate_bn!(BN_mlp)

#     println("########### iteration $(i) completed ###############")
# end

std = BooleanNetwork[]
fit = BooleanNetwork[]
mlp = BooleanNetwork[]

for i = 1:30
    println("########## Starting iteration $(i) #################")

    BN_std, BN_fit, BN_mlp = getBN!(df_bool, 2.0, expect, universe)

    println("Structural acc for STD")
    evaluate_bn!(BN_std)
    # stru_acc = get_structural_accuracy(expect, BN_std.actual_conn, universe)
    dc1 = dynamic_consistency(df_origin_int, BN_std.actual_conn)
    println("Validated dynamic acc $(dc1)")

    println("Structural acc for FIT")
    evaluate_bn!(BN_fit)
    # stru_acc = get_structural_accuracy(expect, BN_fit.actual_conn, universe)
    dc2 = dynamic_consistency(df_origin_int, BN_fit.actual_conn)
    println("Validated dynamic acc $(dc2)")
    
    println("Structural acc for MLP")
    evaluate_bn!(BN_mlp)
    # stru_acc = get_structural_accuracy(expect, BN_mlp.actual_conn, universe)
    dc3 = dynamic_consistency(df_origin_int, BN_mlp.actual_conn)
    println("Validated dynamic acc $(dc3)")

    push!(std, BN_std)
    push!(fit, BN_fit)
    push!(mlp, BN_mlp)

    println("########### iteration $(i) completed ###############")
end

println("Statistics for STD")
stat_bns!(std)
println("Statistics for FIT")
stat_bns!(fit)
println("Statistics for MLP")
stat_bns!(mlp)
