function random_modification(tracker :: CorrelationTracker)
    rndidx(shape) = CartesianIndex((rand(1:s) for s in shape)...)
    return tracker |> size |> rndidx
end
