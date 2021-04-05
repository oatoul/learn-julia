using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate

function f(a::Bool, b::Bool, c::Bool, d::Bool)
    return (a && b) || !c 
end

function data_setup()
    N = 4
    X = reverse.(Iterators.product(fill([false, true],N)...))[:]
    Y = []
    for (a, b, c, d) in X
        append!(Y, f(a, b, c, d))
    end
    X, Y
end

X, Y = data_setup()

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    accuracy = 0.0
    for i in 1:size(X, 1)
        out = process(ind, collect(X[i]))
        if out[1] == Y[i]
            accuracy += 1
        end
    end
    [accuracy / size(X, 1)]
end

cfg = get_config("cfg/boolean.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
run!(e)