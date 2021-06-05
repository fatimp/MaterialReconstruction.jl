module MaterialReconstruction
using CorrelationFunctions.Directional
using CorrelationTrackers
using Statistics: mean, std
using Base.Iterators

include("utilities.jl")
include("modifiers.jl")
include("initialization.jl")
include("furnace.jl")
include("cost.jl")
include("cooldown.jl")
include("annealing.jl")

export
    # Cost functions
    euclid_mean,
    euclid_directional,
    euclid_mean_weighted,
    euclid_directional_weighted,
    čapek_cost,
    # Initialization
    initialize_random,
    # Modifier structures
    AbstractModifier,
    RandomSwapper, RandomFlipper,
    InterfaceSwapper, InterfaceFlipper,
    # Cooldown schedules
    exponential_cooldown,
    aarts_korst_cooldown,
    frost_heineman_cooldown,
    # Annealing
    Furnace, annealing_step

end # module
