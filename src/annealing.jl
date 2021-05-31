function sumcost(tracker1 :: CorrelationTracker,
                 tracker2 :: CorrelationTracker,
                 costfn   :: Function)
    @assert tracked_data(tracker1) == tracked_data(tracker2)
    return mapreduce(data -> costfn(data(tracker1), data(tracker2), data),
                     +, tracked_data(tracker1))
end

"""
    annealing_step(furnace; cooldown = 0.99999, cost = euclid_mean[, modifier])

Perform one step of annealing procedure. `cooldown` is a parameter
defining how fast the furnace will loose temperature. `cost`
determines a function we want to minimize. `modifier` determines which
small modifications are made to the system during the step.

**NB:** The updated state of the annealing process is returned by this
function as a new `Furnace` object. Do not discard it.

See also: [`euclid_mean`](@ref), [`euclid_directional`](@ref),
[`euclid_mean_weighted`](@ref), [`euclid_directional_weighted`](@ref),
[`RandomSwapper`](@ref), [`RandomFlipper`](@ref).
"""
function annealing_step(furnace  :: Furnace;
                        cooldown :: Float64          = 0.99999,
                        cost     :: Function         = euclid_mean,
                        modifier :: AbstractModifier = RandomFlipper(furnace.system))
    # Some statistics
    rejected = false
    accepted = false

    # Compute the cost function which we are trying to minimize
    c1 = sumcost(furnace.system, furnace.target, cost)

    # Do a small modification to our system
    rollback = modify!(modifier)

    # Compute the new value for the cost function
    c2 = sumcost(furnace.system, furnace.target, cost)

    # if c1 < c2 the swap is accepted
    if c2 > c1
        threshold = exp(-(c2 - c1) / furnace.temperature)
        if rand(Float64) > threshold
            rollback!(modifier, rollback)
            @assert c1 ≈ sumcost(furnace.system, furnace.target, cost)
            rejected = true
        end

        accepted = !rejected
    end

    # Construct a new furnace
    return Furnace(furnace.system, furnace.target,
                   furnace.temperature * cooldown,
                   furnace.steps + 1,
                   furnace.rejected + rejected,
                   furnace.accepted + accepted)
end
