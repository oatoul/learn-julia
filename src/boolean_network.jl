export Result, BooleanNetwork, add_dynamic_consistency!, evaluate_bn!, stat_bns!
using Base
import Statistics
import Base: String

Base.@kwdef mutable struct Result
    dynamic_acc::Float64 = 0.0
    struct_acc::Float64 = 0.0
    precision::Float64 = 0.0
    recall::Float64 = 0.0
    TP::Int64 = -1
    FP::Int64 = -1
    FN::Int64 = -1
    TN::Int64 = -1
end

function String(res::Result)
    JSON.json(
        Dict(
            "dynamic_acc"=>res.dynamic_acc,
            "struct_acc"=>res.struct_acc,
            "precision"=>res.precision,
            "recall"=>res.recall,
            "TP"=>res.TP,
            "FP"=>res.FP,
            "FN"=>res.FN,
            "TN"=>res.TN
            )
        )
end

Base.@kwdef mutable struct BooleanNetwork
    actual_conn::Set = Set()
    exptect_conn::Set = Set()
    universe_conn::Set = Set()
    dynamic_consistency::Array{Float64} = Float64[]
    result::Result = Result()
end

function stat_bns!(bns::Array{BooleanNetwork})
    st = Float64[]
    dy = Float64[]
    pr = Float64[]
    re = Float64[]
    TP = Int64[]
    FP = Int64[]
    FN = Int64[]
    TN = Int64[]

    for bn in bns
        push!(st, bn.result.struct_acc)
        push!(dy, bn.result.dynamic_acc)
        push!(pr, bn.result.precision)
        push!(re, bn.result.recall)
        push!(TP, bn.result.TP)
        push!(FP, bn.result.FP)
        push!(FN, bn.result.FN)
        push!(TN, bn.result.TN)
    end
    
    println("###############Statistics######################")
    println("Structural accuracy")
    print_stat!(st)
    println("Dynamic accuracy")
    print_stat!(dy)
    println("Precision")
    print_stat!(pr)
    println("Recall")
    print_stat!(re)
    println("TP")
    print_stat!(TP)
    println("FP")
    print_stat!(FP)
    println("FN")
    print_stat!(FN)
    println("TN")
    print_stat!(TN)
    println("###############################################")

end

function print_stat!(arr)
    println("Max: $(maximum(arr)) Min: $(minimum(arr)) Mean: $(Statistics.mean(arr)) Var: $(Statistics.var(arr))")
end

function add_dynamic_consistency!(bn::BooleanNetwork, c::Float64)
    push!(bn.dynamic_consistency, c)    
end

# function summary(io::IO, ind::CGPInd)
#     print(io, string("Boolean Network Result(", get_active_nodes(ind), ", ",
#                      findall([n.active for n in ind.nodes]), ", ",
#                      ind.outputs, " ,",
#                      ind.fitness, ")"))
# end

function evaluate_bn!(bn::BooleanNetwork)
    expect = bn.exptect_conn
    actual = bn.actual_conn
    universe = bn.universe_conn

    TP = length(intersect(actual, expect))
    FP = length(setdiff(actual, expect))
    FN = length(setdiff(expect, actual))
    expect_negative = setdiff(universe, expect)
    actual_negative = setdiff(universe, actual)
    TN = length(intersect(actual_negative, expect_negative))

    struct_acc = (TP + TN) / (TP + FP + FN + TN)
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    dynamic_acc = Statistics.mean(bn.dynamic_consistency)

    bn.result.TP = TP
    bn.result.FN = FN
    bn.result.FP = FP
    bn.result.TN = TN
    bn.result.struct_acc = struct_acc
    bn.result.precision = precision
    bn.result.recall = recall
    bn.result.dynamic_acc = dynamic_acc

    println("TP : $(TP)")
    println("FP : $(FP)")
    println("FN : $(FN)")
    println("TN : $(TN)")
    println("Structural accuracy : $(struct_acc)")
    println("Dynamic accuracy : $(dynamic_acc)")
    println("Precision : $(precision)")
    println("Recall : $(recall)")
    
    bn
end