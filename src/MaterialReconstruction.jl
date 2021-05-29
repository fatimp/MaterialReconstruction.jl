module MaterialReconstruction
using CorrelationFunctions.Directional
using CorrelationTrackers
using StatsBase: mean

include("modifications.jl")
include("furnace.jl")
include("cost.jl")
include("annealing.jl")

export
    # Cost functions
    euclid_mean,
    euclid_directional,
    euclid_mean_weighted,
    euclid_directional_weighted,
    # Modification functions
    random_modification,
    # Annealing
    Furnace,
    annealing_step

end # module
