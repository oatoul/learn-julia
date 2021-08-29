export Result, BooleanNetwork, add_dynamic_consistency!, evaluate_bn!
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

    bn.result.TP = TP
    bn.result.FN = FN
    bn.result.FP = FP
    bn.result.TN = TN
    bn.result.struct_acc = struct_acc
    bn.result.precision = precision
    bn.result.recall = recall
    bn.result.dynamic_acc = Statistics.mean(bn.dynamic_consistency)

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