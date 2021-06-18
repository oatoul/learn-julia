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

function get_active_connections(ind::CGPInd, low::Int64, high::Int64)
    nodes = ind.nodes[[n.active for n in ind.nodes]]
    res = Set()
    for node in nodes
        x = node.x
        y = node.y
            
        if(x >= low && x <= high)
            push!(res, x)
        end
        if(y >= low && y <= high)
            push!(res, y)
        end
    end
    res
end

df = DataFrame(CSV.File("data/CDC15_bool.tsv",drop=["Time"],type=Bool))
insertcols!(df, 1, :T0 => false, :T1 => true)

ndf = names(df)
l_idx = 3
h_idx = size(ndf)[1]

for col in ndf
    println("$(col)")
end

function runBN!(e::AbstractEvolution, low::Int64, high::Int64, fitness::Int64)
    step!(e)
    while(e.population[end].fitness[1] < fitness && e.gen < e.config.n_gen)
        step!(e)
    end
    log_gen(e)
    save_gen(e)
    set = get_active_connections(e.population[end], low, high)
    println(set)
end

function getBN!(df::DataFrame, fitness::Float64)
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

        
        fit(ind::CGPInd) = evaluate(ind, X, Y)
        e = CGPEvolution(cfg, fit)

        step!(e)
        while(e.population[end].fitness[1] < fitness && e.gen < e.config.n_gen)
            step!(e)
        end

        println("Fitness: $(e.population[end].fitness[1])")
        # println(summary(e.population[end]))

        set = get_active_connections(e.population[end], l_idx, h_idx)
        res = []
        for j in set
            push!(res, ndf[j])
        end

        println("$(target) is connected with:")
        println(res)

    end    
end

cfg = get_config("cfg/CDC15.yaml")
mutate(ind::CGPInd) = goldman_mutate(cfg, ind)
getBN!(df, 0.9)