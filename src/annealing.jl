function sumcost(tracker1 :: CorrelationTracker,
                 tracker2 :: CorrelationTracker,
                 costfn   :: Function)
    @assert tracked_data(tracker1) == tracked_data(tracker2)
    return mapreduce(data -> costfn(data(tracker1), data(tracker2), data),
                     +, tracked_data(tracker1))
end

function annealing_step(furnace      :: Furnace;
                        cooldown     :: Float64 = 0.99999,
                        cost         :: Function = euclid_mean,
                        modification :: Function = random_modification)
    # Some statistics
    rejected = false
    accepted = false

    # Compute the cost function which we are trying to minimize
    c1 = sumcost(furnace.system, furnace.target, cost)

    # Pick a point and swap phase
    idx = modification(furnace.system)
    furnace.system[idx] = 1 - furnace.system[idx]

    # Compute the new value for the cost function
    c2 = sumcost(furnace.system, furnace.target, cost)

    # if c1 < c2 the modification is accepted
    if c2 > c1
        threshold = exp(-(c2 - c1) / furnace.temperature)
        if rand(Float64) > threshold
            # Poor man's reject.
            furnace.system[idx] = 1 - furnace.system[idx]
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
