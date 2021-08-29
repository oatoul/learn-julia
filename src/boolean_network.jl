export Result, BooleanNetwork, add_dynamic_consistency!

using Base

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

# function evaluate(bn::BooleanNetwork)
#     expect = bn.exptect_conn
#     actual = bn.actual_conn
#     universe = bn.universe_conn

#     TP = length(intersect(actual, expect))
#     FP = length(setdiff(actual, expect))
#     FN = length(setdiff(expect, actual))
#     expect_negative = setdiff(universe, expect)
#     actual_negative = setdiff(universe, actual)
#     TN = length(intersect(actual_negative, expect_negative))

#     acc = (TP + TN) / (TP + FP + FN + TN)

#     bn
# end