module CartesianGeneticProgramming

using Cambrian: isless
using Base: Int16
import YAML
using Cambrian
import JSON

include("functions.jl")
include("config.jl")
include("individual.jl")
include("process.jl")
include("mutation.jl")
include("oneplus.jl")
include("crossover.jl")
include("evolution.jl")
include("selection.jl")
include("step.jl")
include("evaluator.jl")
include("mlp_flux.jl")

end
