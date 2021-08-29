export BooleanNetwork, evaluate

mutable struct BooleanNetwork
    dynamic_acc::Float64
    dynamic_consistency::Array{Float64}
    struct_acc::Float64
    precision::Float64
    recall::Float64
    TP::Int64
    FP::Int64
    FN::Int64
    TN::Int64
    actual_conn::Set
    exptect_conn::Set
    universe_conn::Set
    BooleanNetwork() = new()
end

function evaluate(bn::BooleanNetwork)
    expect = bn.exptect_conn
    actual = bn.actual_conn
    universe = bn.universe_conn

    TP = length(intersect(actual, expect))
    FP = length(setdiff(actual, expect))
    FN = length(setdiff(expect, actual))
    expect_negative = setdiff(universe, expect)
    actual_negative = setdiff(universe, actual)
    TN = length(intersect(actual_negative, expect_negative))

    acc = (TP + TN) / (TP + FP + FN + TN)

    bn.TP = TP
    bn.FP = FP
    bn.FN = FN
    bn.TN = TN
    bn.struct_acc = acc

    println("TP : $(bn.TP)")
    println("FP : $(bn.FP)")
    println("FN : $(bn.FN)")
    println("TN : $(bn.TN)")
    println("Structural accuracy : $(bn.struct_acc)")

    bn
end