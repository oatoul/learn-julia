export Node, CGPInd, get_genes, set_genes!, reset!, forward_connections, get_output_trace, get_active_connections!
import Base: copy, String, show, summary, isless, isequal, hash
import Cambrian.print

"default function for nodes, will cause error if used as a function node"
function f_null(args...)::Nothing
    nothing
end

struct Node
    x::Int16
    y::Int16
    f::Function
    active::Bool
end

struct CGPInd <: Cambrian.Individual
    n_in::Int16
    n_out::Int16
    chromosome::Array{Float64}
    genes::Array{Int16}
    outputs::Array{Int16}
    nodes::Array{Node}
    buffer::Array{Float64}
    fitness::Array{Float64}
    sparsity::Int16
    n_active::Int16
end

# "compare graph sparsity when fitness equals"
# function isless(i1::CGPInd, i2::CGPInd)
#     isless(((i1.fitness[1],i2.fitness[1]),(i2.sparsity,i1.sparsity)),(i2.n_active,i1.n_active))
#     # (i1.fitness[1] < i2.fitness[1]) || (i1.sparsity > i2.sparsity) || (i1.n_active > i2.n_active)
#     # (i1.fitness[1] < i2.fitness[1]) || (i1.sparsity > i2.sparsity) || (i1.n_active > i2.n_active)
#     # res = false
#     # if i1.fitness[1] < i2.fitness[1]
#     #     res = true
#     # elseif i2.sparsity > 0 && i1.sparsity > i2.sparsity
#     #     res = true
#     # elseif i2.sparsity == 0 && i1.sparsity > 0
#     #     res =  false
#     # elseif i2.n_active > 0 && i1.n_active > i2.n_active
#     #     res = true
#     # elseif i2.n_active == 0 && i1.n_active > 0
#     #     res = false
#     # end
#     # res
# end

# isless(i1::CGPInd, i2::CGPInd) = isless((i1.fitness[1],i2.fitness[1]),(i2.sparsity,i1.sparsity),(i2.n_active,i1.n_active))

# function isequal(i1::CGPInd, i2::CGPInd)
#     abs(i1.fitness[1] - i2.fitness[1]) < 0.001 && i1.sparsity == i2.sparsity && i1.n_active == i2.n_active
# end

# function hash(x::CGPInd, h::UInt)
#     # key = "$(round(x.fitness[1], digits=3))_$(x.sparsity)_$(x.n_active)"
#     # key = round(x.fitness[1], digits=3) + x.sparsity + x.n_active*100
#     key = (round(x.fitness[1], digits=3), x.sparsity, x.n_active)
#     hash(key, h)
# end

function get_sparsity!(ind::CGPInd)
    "col 1 & 2 are dummy inputs"
    low_idx = 3
    high_idx = Int64(ind.n_in)
    length(get_active_connections!(ind, low_idx, high_idx))
end

function get_active_connections!(nodes_in::Array{Node}, low::Int64, high::Int64)
    nodes = nodes_in[[n.active for n in nodes_in]]
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

function get_active_connections!(ind::CGPInd, low::Int64, high::Int64)
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

function recur_active(active::BitArray, n_in::Integer, ind::Int16, xs::Array{Int16},
                       ys::Array{Int16}, fs::Array{Int16},
                       two_arity::BitArray)::BitArray
    if ind > 0 && ~active[ind]
        active[ind] = true
        active = recur_active(active, n_in, xs[ind], xs, ys, fs, two_arity)
        if two_arity[fs[ind]]
            active = recur_active(active, n_in, ys[ind], xs, ys, fs, two_arity)
        end
    end
    active
end

function find_active(cfg::NamedTuple, genes::Array{Int16},
                     outputs::Array{Int16})::BitArray
    R = cfg.rows
    C = cfg.columns
    n = Int16(cfg.n_in)
    active = falses(R, C)
    xs = genes[:, :, 1] .- n
    ys = genes[:, :, 2] .- n
    fs = genes[:, :, 3]
    for i in eachindex(outputs)
        active = recur_active(active, cfg.n_in, outputs[i] - n, xs, ys, fs,
                              cfg.two_arity)
    end
    active
end

function CGPInd(cfg::NamedTuple, chromosome::Array{Float64}, genes::Array{Int16},
                outputs::Array{Int16})::CGPInd
    R = cfg.rows
    C = cfg.columns
    nodes = Array{Node}(undef, R * C + cfg.n_in)
    for i in 1:cfg.n_in
        nodes[i] = Node(0, 0, f_null, false)
    end
    i = cfg.n_in
    active = find_active(cfg, genes, outputs)
    for y in 1:C
        for x in 1:R
            i += 1
            nodes[i] = Node(genes[x, y, 1], genes[x, y, 2],
                            cfg.functions[genes[x, y, 3]],
                            active[x, y])
        end
    end
    buffer = zeros(R * C + cfg.n_in)
    fitness = -Inf .* ones(cfg.d_fitness)

    sparsity = length(get_active_connections!(nodes, cfg.in_idx, cfg.n_in))
    n_active = size(nodes[[n.active for n in nodes]])[1]
    CGPInd(cfg.n_in, cfg.n_out, chromosome, genes, outputs, nodes, buffer, fitness, sparsity, n_active)
end

function CGPInd(cfg::NamedTuple, chromosome::Array{Float64})::CGPInd
    R = cfg.rows
    C = cfg.columns
    # chromosome: node genes, output genes
    genes = reshape(chromosome[1:(R*C*3)], R, C, 3)
    # TODO: recurrency is ugly and slow
    maxs = collect(1:R:R*C)
    maxs = round.((R*C .- maxs) .* cfg.recur .+ maxs)
    maxs = min.(R*C + cfg.n_in, maxs .+ cfg.n_in)
    maxs = repeat(maxs, 1, R)'
    genes[:, :, 1] .*= maxs
    genes[:, :, 2] .*= maxs
    genes[:, :, 3] .*= length(cfg.functions)
    genes = Int16.(ceil.(genes))
    outputs = Int16.(ceil.(chromosome[(R*C*3+1):end] .* (R * C + cfg.n_in)))
    CGPInd(cfg, chromosome, genes, outputs)
end

function CGPInd(cfg::NamedTuple)::CGPInd
    chromosome = rand(cfg.rows * cfg.columns * 3 + cfg.n_out)
    CGPInd(cfg, chromosome)
end

function CGPInd(cfg::NamedTuple, ind::String)::CGPInd
    dict = JSON.parse(ind)
    CGPInd(cfg, Array{Float64}(dict["chromosome"]))
end

function copy(n::Node)
    Node(n.x, n.y, n.f, n.active)
end

function copy(ind::CGPInd)
    buffer = zeros(length(ind.buffer))
    nodes = Array{Node}(undef, length(ind.nodes))
    for i in eachindex(ind.nodes)
        nodes[i] = copy(ind.nodes[i])
    end
    CGPInd(ind.n_in, ind.n_out, copy(ind.chromosome), copy(ind.genes),
           copy(ind.outputs), nodes, buffer, copy(ind.fitness))
end

function String(n::Node)
    JSON.json(Dict(:x=>n.x, :y=>n.y, :f=>string(n.f), :active=>n.active))
end

function show(io::IO, n::Node)
    print(io, String(n))
end

function String(ind::CGPInd)
    JSON.json(Dict("chromosome"=>ind.chromosome, "fitness"=>ind.fitness,
                   "xs"=>ind.genes[:, :, 1][:], "ys"=>ind.genes[:, :, 2][:],
                   "fs"=>ind.genes[:, :, 3][:]))
end

function show(io::IO, ind::CGPInd)
    print(io, String(ind))
end

function get_active_nodes(ind::CGPInd)
    ind.nodes[[n.active for n in ind.nodes]]
end

function get_active_number(ind::CGPInd)
    size(get_active_nodes(ind))
end

function summary(io::IO, ind::CGPInd)
    print(io, string("CGPInd(", get_active_nodes(ind), ", ",
                     findall([n.active for n in ind.nodes]), ", ",
                     ind.outputs, " ,",
                     ind.fitness, ")"))
end

function interpret(i::CGPInd)
    x::AbstractArray->process(i, x)
end

function reset!(c::CGPInd)
    c.buffer .= 0.0
end

function get_genes(c::CGPInd, node_id::Integer)::Array{Float64}
    if node_id > c.n_in
        return c.chromosome[(node_id-c.n_in-1)*3 .+ (1:3)]
    else
        return zeros(3)
    end
end

function get_genes(c::CGPInd, nodes::Array{<:Integer})::Array{Float64}
    if length(nodes) > 0
        return reduce(vcat, map(x->get_genes(c, x), nodes))
    else
        return Array{Float64}(0)
    end
end

"set the genes of node_id to genes"
function set_genes!(c::CGPInd, node_id::Integer, genes::Array{Float64})
    if node_id > c.n_in
        @assert length(genes) == 3
        c.chromosome[(node_id-c.n_in-1)*3 .+ (1:3)] = genes
    end
end

"a list for each node i of a list of all the nodes which have the node i as an input, and i"
function forward_connections(c::CGPInd)
    connections = [[i] for i in 1:length(c.nodes)]
    for ci in eachindex(c.nodes)
        conns = [c.nodes[ci].x, c.nodes[ci].y]
        for conn in conns
            if conn > 0
                push!(connections[conn], ci)
            end
        end
    end
    map(unique, connections)
end

"recursively walk back through inputs"
function recur_output_trace(c::CGPInd, ind::Int16, visited::Array{Int16})
    if ~(ind in visited)
        push!(visited, ind)
        if ind > c.n_in
            recur_output_trace(c, c.nodes[ind].x, visited)
            if CGPFunctions.arity[string(c.nodes[ind].f)] == 2
                recur_output_trace(c, c.nodes[ind].y, visited)
            end
        end
    end
    visited
end

"return a list of node indices which determine the output"
function get_output_trace(c::CGPInd, output_ind::Integer)
    recur_output_trace(c, c.outputs[output_ind], Array{Int16}(undef, 0))
end

"return a list of node indices which determine the output for multiple outputs"
function get_output_trace(c::CGPInd, outputs::Array{<:Integer})
    if length(outputs) > 0
        return unique(reduce(vcat, map(x->get_output_trace(c, x), outputs)))
    else
        return Array{Int16}(undef, 0)
    end
end

"return a list of node indices which determine the output for all outputs"
function get_output_trace(c::CGPInd)
    get_output_trace(c, collect(1:c.n_out))
end
