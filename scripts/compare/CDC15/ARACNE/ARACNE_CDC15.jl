using CSV
using DataFrames
using CartesianGeneticProgramming

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

function get_actual()
    res = Set()
    push!(res, "ACE2_SWI5")
    push!(res, "MBP1_NDD1")
    push!(res, "NDD1_MBP1")
    push!(res, "SWI6_MBP1")
    push!(res, "SWI4_SWI5")
    push!(res, "SWI5_MBP1")
    push!(res, "ACE2_NDD1")
    push!(res, "SWI5_ACE2")
    push!(res, "SWI5_SWI4")
    push!(res, "MBP1_SWI6")
    push!(res, "MBP1_SWI5")
    push!(res, "MBP1_SWI4")
    push!(res, "SWI4_MBP1")
    push!(res, "NDD1_ACE2")
    
    res
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

df_origin = DataFrame(CSV.File("data/CDC15_bool.tsv",drop=["Time"],type=Int))

universe = get_universe_set(names(df_origin))

expect = get_expect_CDC15()

actual = get_actual()

stru_acc = get_structural_accuracy(expect, actual, universe)

dc = dynamic_consistency(df_origin, get_actual())

println("Calculated dynamic accuracy: $(dc)")