metric(x :: AbstractArray, y :: AbstractArray) =
    mapreduce((x, y) -> (x - y)^2, +, x, y)

# All cost functions accept three arguments: the first two of type
# CorrelationData and the last of type TrackedData (which holds the
# information about what is measured in the first two args and can be
# ignored).

euclid_mean(data1 :: CorrelationData, data2 :: CorrelationData) =
    metric(data1 |> mean, data2 |> mean)

euclid_mean(data1 :: CorrelationData, data2 :: CorrelationData, :: TrackedData) =
    euclid_mean(data1, data2)

function euclid_directional(data1 :: CorrelationData,
                            data2 :: CorrelationData)
    @assert directions(data1) == directions(data2)
    sum(metric(data1[dir], data2[dir]) for dir in directions(data1))
end

euclid_directional(data1 :: CorrelationData, data2 :: CorrelationData, :: TrackedData) =
    euclid_directional(data1, data2)

function euclid_mean_weighted(furnace :: Furnace)
    @assert tracked_data(furnace.system) == tracked_data(furnace.target)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_mean(furnace.system |> data,
                                       furnace.target |> data)
                   for data in tracked_data(furnace.target))
    return (data1    :: CorrelationData,
            data2    :: CorrelationData,
            datatype :: TrackedData) -> euclid_mean(data1, data2) / weights[datatype]
end

function euclid_directional_weighted(furnace :: Furnace)
    @assert tracked_data(furnace.system) == tracked_data(furnace.target)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_directional(furnace.system |> data,
                                              furnace.target |> data)
                   for data in tracked_data(furnace.target))
    return (data1    :: CorrelationData,
            data2    :: CorrelationData,
            datatype :: TrackedData) -> euclid_directional(data1, data2) / weights[datatype]
end
