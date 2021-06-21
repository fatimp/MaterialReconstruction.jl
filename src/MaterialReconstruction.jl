module MaterialReconstruction
using CorrelationFunctions.Directional
using CorrelationTrackers
using Base.Iterators
using Statistics: mean, std
using LsqFit: curve_fit, coef
using PoissonRandom: pois_rand
using CircularArrays: CircularArray
using ImageSegmentation: label_components

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
    initialize_spheres,
    # Modifier structures
    AbstractModifier,
    RandomSwapper, RandomFlipper,
    InterfaceSwapper, InterfaceFlipper,
    C2UpdateProxy,
    modify!, rollback!,
    # Cooldown schedules
    exponential_cooldown,
    aarts_korst_cooldown,
    frost_heineman_cooldown,
    # Annealing
    Furnace, annealing_step

end # module
