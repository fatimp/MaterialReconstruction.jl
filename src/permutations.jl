function random_permutation(tracker :: CorrelationTracker)
    shape = size(tracker)
    rndidx(shape) = CartesianIndex((rand(1:s) for s in shape)...)

    idx1 = rndidx(shape)
    while true
        # Try to find element with a different phase
        idx2 = rndidx(shape)
        if tracker[idx1] != tracker[idx2]
            return idx1, idx2
        end
    end
end
