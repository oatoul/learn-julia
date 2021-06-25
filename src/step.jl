export step!

import Cambrian.step!
import Cambrian.AbstractEvolution

"""
update e.elites with the best individuals from the current population or
existing elites, based on fitness
"""
function elites_generation(e::AbstractEvolution)
    pop = sort(unique([e.population; e.elites]), by=i ->(i.fitness, -i.sparsity, -i.n_active))
    # for el in pop
    #     println("Fitness $(el.fitness) sparsity $(el.sparsity) active $(el.n_active)")
    # end
    # println("...................................................")
    n_elites = length(e.elites)
    e.elites = deepcopy(pop[end-n_elites+1:end])
end

"""
generation should be used for any additional record keeping necessary in an
algorithm. It is called after populate and evaluate, so all individuals should
have genes and fitness reflecting their current state. The default function
calls elites_generation if e.elites is defined, otherwise does nothing.
"""
function generation(e::AbstractEvolution)
    if hasfield(typeof(e), :elites)
        elites_generation(e)
    end
end

"""
The generic iteration of an evolution. Calls populate, evaluate, and generation.
Also calls log_gen and save_gen based on the provided config values. Subclasses
of AbstractEvolution should override the populate, evaluate, or generation
functions rather than overriding this function.
"""
function step!(e::AbstractEvolution)
    e.gen += 1
    if e.gen > 1
        populate(e)
    end
    evaluate(e)
    generation(e)

    if ((e.config.log_gen > 0) && mod(e.gen, e.config.log_gen) == 0)
        log_gen(e)
    end
    if ((e.config.save_gen > 0) && mod(e.gen, e.config.save_gen) == 0)
        save_gen(e)
    end
end

"Call step!(e) e.config.n_gen times consecutively"
function run!(e::AbstractEvolution)
    for i in (e.gen+1):e.config.n_gen
        step!(e)
    end
end
