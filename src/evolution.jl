export CGPEvolution, save_gen

import Cambrian.evaluate

mutable struct CGPEvolution{T} <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{T}
    fitness::Function
    gen::Int
end

populate(e::CGPEvolution) = oneplus_populate_v2(e)
evaluate(e::CGPEvolution) = Cambrian.fitness_evaluate(e, e.fitness)

function CGPEvolution(cfg::NamedTuple, fitness::Function;
                      logfile=string("logs/", replace(cfg.id,":" => "_"), ".csv"))
    logger = CambrianLogger(logfile)
    population = Cambrian.initialize(CGPInd, cfg)
    CGPEvolution(cfg, logger, population, fitness, 0)
end
