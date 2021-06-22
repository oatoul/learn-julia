using CartesianGeneticProgramming
using Cambrian
using CSV
using DataFrames
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
        e = CGPEvolution(cfg, fit)

        step!(e)
        while(e.population[end].fitness[1] < fitness && e.gen < e.config.n_gen)
            step!(e)
        end

        println("Fitness: $(e.population[end].fitness[1])")
        # println(summary(e.population[end]))

        set = get_active_connections!(e.population[end], l_idx, h_idx)

        "output result for target gene"
        # res = []
        # for j in set
        #     push!(res, ndf[j])
        # end
        # println("$(target) is connected with:")
        # println(res)

        "store connections to BN"
        for k in set
            conn = ndf[k] * "_" * target
            push!(BN, conn)
            println(conn)
        end


    end

    BN
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

# ndf = names(df)
# l_idx = 3
# h_idx = size(ndf)[1]
# for col in ndf
#     println("$(col)")
# end

universe = get_universe_set(names(df_origin))

expect = get_expect_CDC15()

actual = getBN!(df, 0.95)

stru_acc = get_structural_accuracy(expect, actual, universe)
