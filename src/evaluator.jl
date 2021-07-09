export ConfusionMatrix, compute_confusion_matrix
# Alias for the NamedTuple type
const ConfusionMatrix = NamedTuple{(:TP, :FP, :TN, :FN),NTuple{4,Int64}}

"""
    compute_confusion_matrix(true_net::AbstractMatrix{Bool}, 
        inferred_net::AbstractMatrix{Bool})::ConfusionMatrix

The (i, j)-th element in a network matrix denotes whether a link exists between the i-th node and 
the j-th node. By default, the links are `directed`.
"""
function compute_confusion_matrix(true_net::AbstractMatrix{Bool}, 
    inferred_net::AbstractMatrix{Bool}; directed::Bool=true)::ConfusionMatrix
    @assert size(true_net) == size(inferred_net)
    tp = fp = tn = fn = 0
    for (x, y) in zip(true_net, inferred_net)
        if x
            if y
                tp += 1
            else
                fn += 1
            end
        else
            if y 
                fp += 1
            else
                tn += 1
            end
        end
    end
    # if the edge is undirected, then each wrong/right case is counted twice.
    if !directed
        tp = tp ÷ 2
        fp = fp ÷ 2
        tn = tn ÷ 2
        fn = fn ÷ 2
    end
    return (TP = tp, FP = fp, TN = tn, FN = fn)
end

# example
# cm = compute_confusion_matrix([true false; false true], [false true; true true])
# println(cm.FP, ", ", cm.FN)
