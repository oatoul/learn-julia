export CGPEvolution, save_gen

import Cambrian.evaluate

mutable struct CGPEvolution{T} <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{T}
    fitness::Function
    gen::Int
    elites::Array{T}
end

populate(e::CGPEvolution) = oneplus_populate_v2(e)
evaluate(e::CGPEvolution) = Cambrian.fitness_evaluate(e, e.fitness)

function CGPEvolution(cfg::NamedTuple, fitness::Function;
                      logfile=string("logs/", replace(cfg.id,":" => "_"), ".csv"))
    logger = CambrianLogger(logfile)
    population = Cambrian.initialize(CGPInd, cfg)
    elites = initialize_elites(CGPInd, cfg)
    # println("Initialize CGPEvolution")
    CGPEvolution(cfg, logger, population, fitness, 0, elites)
end

"create all members of the first elites"
function initialize_elites(itype::Type, cfg::NamedTuple)
    # println("Initialize n_elite = $(cfg.n_elite)")
    population = Array{itype}(undef, cfg.n_elite)
    for i in 1:cfg.n_elite
        population[i] = itype(cfg)
    end
    population
end
