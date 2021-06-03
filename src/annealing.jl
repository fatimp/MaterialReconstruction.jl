"""
    annealing_step(furnace[; cooldown][, cost][, modifier])

Perform one step of annealing procedure. `cooldown` is a function
defining how fast the furnace will loose temperature. `cost`
determines a function we want to minimize. `modifier` determines which
small modifications are made to the system during the step.

**NB:** The updated state of the annealing process is returned by this
function as a new `Furnace` object. Do not discard it.

See also: [`euclid_mean`](@ref), [`euclid_directional`](@ref),
[`euclid_mean_weighted`](@ref), [`euclid_directional_weighted`](@ref),
[`RandomSwapper`](@ref), [`RandomFlipper`](@ref),
[`InterfaceSwapper`](@ref), [`InterfaceFlipper`](@ref),
[`exponential_cooldown`](@ref).
"""
function annealing_step(furnace  :: Furnace;
                        cooldown :: Function         = exponential_cooldown(),
                        cost     :: Function         = euclid_mean,
                        modifier :: AbstractModifier = RandomFlipper(furnace.system))
    # Some statistics
    rejected = false
    accepted = false

    # Compute the cost function which we are trying to minimize
    c1 = cost(furnace.system, furnace.target)

    # Do a small modification to our system
    rollback = modify!(modifier)

    # Compute the new value for the cost function
    c2 = cost(furnace.system, furnace.target)

    # if c1 < c2 the swap is accepted
    if c2 > c1
        threshold = exp(-(c2 - c1) / furnace.temperature)
        if rand(Float64) > threshold
            rollback!(modifier, rollback)
            @assert c1 ≈ cost(furnace.system, furnace.target)
            rejected = true
        end

        accepted = !rejected
    end

    # Calculate a new temperature
    newtemperature = rejected ? furnace.temperature : cooldown(furnace.temperature, c2)

    # Construct a new furnace
    return Furnace(furnace.system, furnace.target,
                   newtemperature,
                   furnace.steps + 1,
                   furnace.rejected + rejected,
                   furnace.accepted + accepted)
end
