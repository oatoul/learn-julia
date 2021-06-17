using DataFrames: insert_single_column!
using CartesianGeneticProgramming
using Cambrian
using CSV
using DataFrames
import Cambrian.mutate


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

df = DataFrame(CSV.File("data/CDC15_bool.tsv",drop=["Time"],type=Bool))
insertcols!(df, 1, :T0 => false, :T1 => true)

for col in names(df)
    println("$(col)")
end

function runBN!(e::AbstractEvolution)
    step!(e)
    while(e.population[end].fitness[1] < 0.9 && e.gen < e.config.n_gen)
        step!(e)
    end
    println(summary(e.population[end]))
    log_gen(e)
    save_gen(e)
end

target = "ACE2"

X = copy(Tuple.(eachrow(df)))
Y = copy(df[!,target])
pop!(X)
popfirst!(Y)

cfg = get_config("cfg/CDC15.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
runBN!(e)