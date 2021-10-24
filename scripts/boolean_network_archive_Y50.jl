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


# function runBN!(e::AbstractEvolution, low::Int64, high::Int64, fitness::Int64)
#     step!(e)
#     while(e.population[end].fitness[1] < fitness && e.gen < e.config.n_gen)
#         step!(e)
#     end
#     log_gen(e)
#     save_gen(e)
#     set = get_active_connections(e.population[end], low, high)
#     println(set)
# end

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
        # println("Calculate $(i): $(target)")

        "data setup for target gene"
        Y = copy(df[!,target])
        popfirst!(Y)

        
        fit(ind::CGPInd) = evaluate!(ind, X, Y)
        e = CGPEvolution_archive(cfg, fit)

        step_archive!(e)
        while(e.gen < e.config.n_gen)
            # println("Step Fitness $(e.population[end].fitness[1]) sparsity $(e.population[end].sparsity)")
            step_archive!(e)
        end

        # println("Fitness: $(e.elites[end].fitness[1])")

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

    # println(BN_std.dynamic_consistency)
    # println(BN_fit.dynamic_consistency)
    # println(BN_mlp.dynamic_consistency)

    BN_std, BN_fit, BN_mlp
end


# function store_BN(set::Set, BN::Set, ndf::Array, target::String)
#     for k in set
#         conn = ndf[k] * "_" * target
#         push!(BN, conn)
#     end
# end


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
    # println(res_all)
    count = StatsBase.countmap(res_all)
    # println(count)
    for (k,v) in count
        if v >= min_vote
            push!(res_set,k)
        end
    end
    # println(res_set)
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

    # println(res_all)

    best = sort(res_all, by=((x,y),)->(-x,))[end]
    loss = best[1]
    best = best[2]
    # println("Best: Fitness $(best.fitness) sparsity $(best.sparsity) active $(best.n_active) loss $(loss)")

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

function get_expect()
    res = Set()
    push!(res, "YCL066W_YKL178C")
    push!(res, "YCL067C_YCL066W")
    push!(res, "YCL067C_YGR044C")
    push!(res, "YCR040W_YKL178C")
    push!(res, "YCR096C_YCR039C")
    push!(res, "YCR096C_YCR040W")
    push!(res, "YCR096C_YKL178C")
    push!(res, "YCR097W_YCL066W")
    push!(res, "YCR097W_YCR039C")
    push!(res, "YCR097W_YCR040W")
    push!(res, "YCR097W_YGR044C")
    push!(res, "YCR097W_YKL178C")
    push!(res, "YDR146C_YBL111C")
    push!(res, "YDR146C_YCR018C")
    push!(res, "YDR146C_YDR543C")
    push!(res, "YDR146C_YDR545W")
    push!(res, "YDR146C_YER188CA")
    push!(res, "YDR146C_YER189W")
    push!(res, "YDR146C_YER190W")
    push!(res, "YDR146C_YGR044C")
    push!(res, "YDR146C_YGR296W")
    push!(res, "YDR146C_YLR013W")
    push!(res, "YDR146C_YMR135WA")
    push!(res, "YDR146C_YNL336W")
    push!(res, "YDR146C_YOR140W")
    push!(res, "YDR146C_YPL283C")
    push!(res, "YDR421W_YCL065W")
    push!(res, "YDR421W_YCL066W")
    push!(res, "YDR421W_YCL067C")
    push!(res, "YGL035C_YBR019C")
    push!(res, "YGL035C_YBR020W")
    push!(res, "YGL035C_YDR009W")
    push!(res, "YGL035C_YDR146C")
    push!(res, "YGL035C_YKL109W")
    push!(res, "YGL035C_YMR280C")
    push!(res, "YGL035C_YPL248C")
    push!(res, "YGL209W_YBR019C")
    push!(res, "YGL209W_YBR020W")
    push!(res, "YGL209W_YDR009W")
    push!(res, "YGL209W_YDR146C")
    push!(res, "YGL209W_YGL035C")
    push!(res, "YGL209W_YKL109W")
    push!(res, "YGL209W_YMR280C")
    push!(res, "YGL209W_YPL248C")
    push!(res, "YGR044C_YCR018C")
    push!(res, "YGR044C_YOR140W")
    push!(res, "YHR006W_YCR039C")
    push!(res, "YHR006W_YCR040W")
    push!(res, "YHR006W_YCR041W")
    push!(res, "YJL056C_YCL065W")
    push!(res, "YJL056C_YCL066W")
    push!(res, "YJL056C_YCL067C")
    push!(res, "YJL056C_YCR039C")
    push!(res, "YJL056C_YCR040W")
    push!(res, "YJL056C_YCR041W")
    push!(res, "YJL056C_YDR545W")
    push!(res, "YJL056C_YER188CA")
    push!(res, "YJL056C_YER189W")
    push!(res, "YJL056C_YNL336W")
    push!(res, "YJL056C_YNL337W")
    push!(res, "YJL056C_YNL339C")
    push!(res, "YJL110C_YBL111C")
    push!(res, "YJL110C_YCL065W")
    push!(res, "YJL110C_YCL066W")
    push!(res, "YJL110C_YCL067C")
    push!(res, "YJL110C_YCR039C")
    push!(res, "YJL110C_YCR040W")
    push!(res, "YJL110C_YCR041W")
    push!(res, "YJL110C_YCR096C")
    push!(res, "YJL110C_YDR543C")
    push!(res, "YJL110C_YDR544C")
    push!(res, "YJL110C_YDR545W")
    push!(res, "YJL110C_YER189W")
    push!(res, "YJL110C_YER190W")
    push!(res, "YJL110C_YGL001C")
    push!(res, "YJL110C_YGR296W")
    push!(res, "YJL110C_YHR091C")
    push!(res, "YJL110C_YLR463C")
    push!(res, "YJL110C_YLR465C")
    push!(res, "YJL110C_YLR467W")
    push!(res, "YJL110C_YNL337W")
    push!(res, "YJL110C_YNL339C")
    push!(res, "YKL109W_YBL111C")
    push!(res, "YKL109W_YCL065W")
    push!(res, "YKL109W_YCL066W")
    push!(res, "YKL109W_YCL067C")
    push!(res, "YKL109W_YCR039C")
    push!(res, "YKL109W_YCR040W")
    push!(res, "YKL109W_YCR041W")
    push!(res, "YKL109W_YDR543C")
    push!(res, "YKL109W_YDR544C")
    push!(res, "YKL109W_YDR545W")
    push!(res, "YKL109W_YER188CA")
    push!(res, "YKL109W_YER189W")
    push!(res, "YKL109W_YGL001C")
    push!(res, "YKL109W_YGR296W")
    push!(res, "YKL109W_YHR091C")
    push!(res, "YKL109W_YLR463C")
    push!(res, "YKL109W_YLR465C")
    push!(res, "YKL109W_YLR467W")
    push!(res, "YKL109W_YNL336W")
    push!(res, "YKL109W_YNL337W")
    push!(res, "YKL109W_YNL339C")
    push!(res, "YKL109W_YPL283C")
    push!(res, "YLR013W_YBL111C")
    push!(res, "YLR013W_YCL065W")
    push!(res, "YLR013W_YCL066W")
    push!(res, "YLR013W_YCL067C")
    push!(res, "YLR013W_YCR039C")
    push!(res, "YLR013W_YCR040W")
    push!(res, "YLR013W_YCR041W")
    push!(res, "YLR013W_YDR543C")
    push!(res, "YLR013W_YDR544C")
    push!(res, "YLR013W_YDR545W")
    push!(res, "YLR013W_YER188CA")
    push!(res, "YLR013W_YER189W")
    push!(res, "YLR013W_YGR296W")
    push!(res, "YLR013W_YHR091C")
    push!(res, "YLR013W_YLR463C")
    push!(res, "YLR013W_YLR465C")
    push!(res, "YLR013W_YLR467W")
    push!(res, "YLR013W_YNL336W")
    push!(res, "YLR013W_YNL337W")
    push!(res, "YLR013W_YNL339C")
    push!(res, "YLR013W_YPL283C")
    push!(res, "YLR098C_YCL065W")
    push!(res, "YLR098C_YCL066W")
    push!(res, "YLR098C_YCL067C")
    push!(res, "YLR098C_YCR039C")
    push!(res, "YLR098C_YCR040W")
    push!(res, "YLR098C_YCR041W")
    push!(res, "YLR098C_YER189W")
    push!(res, "YLR098C_YER190W")
    push!(res, "YLR098C_YMR135WA")
    push!(res, "YML051W_YBR019C")
    push!(res, "YML051W_YBR020W")
    push!(res, "YML113W_YCR096C")
    push!(res, "YML113W_YDR543C")
    push!(res, "YML113W_YDR544C")
    push!(res, "YML113W_YDR545W")
    push!(res, "YML113W_YER189W")
    push!(res, "YML113W_YER190W")
    push!(res, "YML113W_YGL001C")
    push!(res, "YML113W_YGR296W")
    push!(res, "YML113W_YHR091C")
    push!(res, "YML113W_YLR463C")
    push!(res, "YML113W_YLR465C")
    push!(res, "YML113W_YLR467W")
    push!(res, "YML113W_YNL337W")
    push!(res, "YML113W_YNL339C")
    push!(res, "YML113W_YPL283C")
    push!(res, "YOL089C_YCR096C")
    push!(res, "YOL089C_YCR097W")
    push!(res, "YOL089C_YGL001C")
    push!(res, "YOL089C_YNL336W")
    push!(res, "YOL089C_YNL337W")
    push!(res, "YOL089C_YNL339C")
    push!(res, "YOR032C_YGL001C")
    push!(res, "YOR032C_YHR006W")
    push!(res, "YPL248C_YBL111C")
    push!(res, "YPL248C_YBR019C")
    push!(res, "YPL248C_YBR020W")
    push!(res, "YPL248C_YCL065W")
    push!(res, "YPL248C_YCL066W")
    push!(res, "YPL248C_YCL067C")
    push!(res, "YPL248C_YDR009W")
    push!(res, "YPL248C_YDR544C")
    push!(res, "YPL248C_YDR545W")
    push!(res, "YPL248C_YER188CA")
    push!(res, "YPL248C_YER189W")
    push!(res, "YPL248C_YER190W")
    push!(res, "YPL248C_YGR296W")
    push!(res, "YPL248C_YHR091C")
    push!(res, "YPL248C_YML051W")
    push!(res, "YPL248C_YMR135WA")
    push!(res, "YPL248C_YNL336W")
    push!(res, "YPL248C_YNL337W")
    push!(res, "YPL248C_YNL339C")
    push!(res, "YPL248C_YOR140W")
    push!(res, "YPR199C_YCR039C")
    push!(res, "YPR199C_YCR040W")
    push!(res, "YPR199C_YCR041W")
    push!(res, "YPR199C_YCR096C")
    push!(res, "YPR199C_YMR135WA")    
    res
end


cfg = get_config("cfg/Yeast50.yaml")
mutate(ind::CGPInd) = goldman_mutate(cfg, ind)

print(cfg)

df_origin_bool = DataFrame(CSV.File("data/Yeast50_bool.tsv",drop=["Time"],type=Bool))
df_origin_int = DataFrame(CSV.File("data/Yeast50_bool.tsv",drop=["Time"],type=Int))
df_bool = copy(df_origin_bool)
insertcols!(df_bool, 1, :T0 => false, :T1 => true)

println(df_bool)

df_origin_float = DataFrame(CSV.File("data/Yeast50.tsv",drop=["Time"],type=Float32))
df_float = Matrix(df_origin_float)'

# ndf = names(df)
# l_idx = 3
# h_idx = size(ndf)[1]
# for col in ndf
#     println("$(col)")
# end

universe = get_universe_set(names(df_origin_bool))

expect = get_expect()

# df_mean = get_mean("data/CDC15.tsv","Time")

# println("Mean is $(df_mean)")

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

# BN_std, BN_fit, BN_mlp = getBN!(df, 0.99)

# # actual = getBN!(df, 0.99)
# println("Structural acc for STD")
# stru_acc = get_structural_accuracy(expect, BN_std, universe)

# println("Structural acc for FIT")
# stru_acc = get_structural_accuracy(expect, BN_fit, universe)

# println("Structural acc for MLP")
# stru_acc = get_structural_accuracy(expect, BN_mlp, universe)
