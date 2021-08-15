function archive_insert(ind::CGPInd, pop::Array{<:CGPInd})
    res = CGPInd[]
    survive = false

    if ind.sparsity == 0
        exist = false
        for p in pop
            if p.sparsity == 0
                exist = true
                if p.fitness > ind.fitness
                    push!(res, p)
                else
                    survive = true
                end
            else
                push!(res, p)
            end
        end

        if exist == false
            survive = true
        end
    else
        for p in pop
            if p.fitness > ind.fitness || p.sparsity < ind.sparsity
                push!(res, p)
            else
                survive = true
            end
        end
    end

    if survive
        push!(res, ind)
    end

    res
end