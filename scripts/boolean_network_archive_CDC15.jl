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

function getBN!(df::DataFrame, fitness::Float64)
    BN = Set()

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
            # println("Step Fitness $(e.population[end].fitness[1]) sparsity $(e.population[end].sparsity)")
            step_archive!(e)
        end

        println("Fitness: $(e.elites[end].fitness[1])")
        # println(summary(e.population[end]))

        # eee = sort(e.elites)
        for el in e.elites
            println("Fitness $(el.fitness) sparsity $(el.sparsity) active $(el.n_active)")
        end

        if size(e.elites)[1] > 1
            set = get_connections_from_archive!(e.elites, l_idx, h_idx, i)
        else
            set = get_active_connections!(e.elites[end], l_idx, h_idx)
        end

        "output result for target gene"
        # res = []
        # for j in set
        #     push!(res, ndf[j])
        # end
        # println("$(target) is connected with:")
        # println(res)

        "store connections to BN"
        # for k in set
        #     conn = ndf[k] * "_" * target
        #     push!(BN, conn)
        #     println(conn)
        # end
        store_BN(set, BN, ndf, target)

    end

    BN
end


function store_BN(set::Set, BN::Set, ndf::Array, target::String)
    for k in set
        conn = ndf[k] * "_" * target
        push!(BN, conn)
        # println(conn)
    end
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

function get_connections_from_archive!(elites::Array{CGPInd}, low::Int64, high::Int64, target_idx::Int64)
    res_all = []
    shift = low - 1
    for e in elites
        err = evaluate_mlp(e, low, high, target_idx-shift)
        res_all = vcat(res_all, (err, e))
        # println(e)
        println("Err: $(err)")
    end

    # println(res_all)

    best = sort(res_all, by=((x,y),)->(-x,))[end]
    err = best[1]
    best = best[2]
    println("Best: Fitness $(best.fitness) sparsity $(best.sparsity) active $(best.n_active) err $(err)")

    get_active_connections!(best, low, high)
end

function evaluate_mlp(ind::CGPInd, low::Int64, high::Int64, target_idx::Int64)
    in_idx = collect(Int64, get_active_connections_no_shift!(ind, low, high))
    err = train_mlp_flux(in_idx, target_idx)
    err
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

function get_expect_CDC15()
    res = Set()
    push!(res, "FKH2_SWI5")
    push!(res, "FKH2_ACE2")

    push!(res, "ACE2_CLN3")

    push!(res, "SWI6_NDD1")
    push!(res, "SWI6_SWI4")

    push!(res, "NDD1_ACE2")
    push!(res, "NDD1_SWI5")

    push!(res, "SWI5_CLN3")

    push!(res, "MCM1_ACE2")
    push!(res, "MCM1_SWI5")
    push!(res, "MCM1_CLN3")
    push!(res, "MCM1_SWI4")

    push!(res, "CLN3_SWI4")
    push!(res, "CLN3_MBP1")

    push!(res, "MBP1_NDD1")
    push!(res, "MBP1_SWI4")

    push!(res, "SWI4_NDD1")

    res
end


cfg = get_config("cfg/CDC15.yaml")
mutate(ind::CGPInd) = goldman_mutate(cfg, ind)

df_origin = DataFrame(CSV.File("data/CDC15_bool.tsv",drop=["Time"],type=Bool))
df = copy(df_origin)
insertcols!(df, 1, :T0 => false, :T1 => true)

println(df)

# ndf = names(df)
# l_idx = 3
# h_idx = size(ndf)[1]
# for col in ndf
#     println("$(col)")
# end

universe = get_universe_set(names(df_origin))

expect = get_expect_CDC15()

BN = getBN!(df, 0.99)

# actual = getBN!(df, 0.99)
println("Structural acc")
stru_acc = get_structural_accuracy(expect, BN, universe)
