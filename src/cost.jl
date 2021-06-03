metric(x :: AbstractArray, y :: AbstractArray) =
    mapreduce((x, y) -> (x - y)^2, +, x, y)

function sumcosts(tracker1 :: CorrelationTracker,
                  tracker2 :: CorrelationTracker,
                  costfn   :: Function,
                           :: Nothing)
    @assert tracked_data(tracker1) == tracked_data(tracker2)
    return mapreduce(data -> costfn(tracker1 |> data,
                                    tracker2 |> data),
                     +, tracked_data(tracker1))
end

function sumcosts(tracker1 :: CorrelationTracker,
                  tracker2 :: CorrelationTracker,
                  costfn   :: Function,
                  weights  :: Dict{TrackedData{N}, Float64}) where N
    @assert tracked_data(tracker1) == tracked_data(tracker2)
    return mapreduce(data -> costfn(tracker1 |> data,
                                    tracker2 |> data) / weights[data],
                     +, tracked_data(tracker1))
end

# All cost functions accept two systems as their arguments: and return
# a distance between those systems. We will then try to minimize that
# distance.

"""
    euclid_mean(data1 :: CorrelationData, data2 :: CorrelationData)
    euclid_mean(data1 :: CorrelationTracker, data2 :: CorrelationTracker)

Calculate euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `CorrelationTracker` object). The values are averaged
along all directions before calculation.
"""
function euclid_mean end

euclid_mean(data1 :: CorrelationData, data2 :: CorrelationData) =
    metric(data1 |> mean, data2 |> mean)

euclid_mean(data1 :: CorrelationTracker, data2 :: CorrelationTracker) =
    sumcosts(data1, data2, euclid_mean, nothing)

"""
    euclid_directional(data1 :: CorrelationData, data2 :: CorrelationData)
    euclid_directional(data1 :: CorrelationTracker, data2 :: CorrelationTracker)

Calculate euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `CorrelationTracker` object). The values calculated
along different directions are not averaged and treated separately.
"""
function euclid_directional end

function euclid_directional(data1 :: CorrelationData,
                            data2 :: CorrelationData)
    @assert directions(data1) == directions(data2)
    sum(metric(data1[dir], data2[dir]) for dir in directions(data1))
end

function euclid_directional(data1 :: CorrelationTracker,
                            data2 :: CorrelationTracker)
    sumcosts(data1, data2, euclid_directional, nothing)
end

"""
    euclid_mean_weighted(data1 :: CorrelationTracker, data2 :: CorrelationTracker)

Calculate euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `CorrelationTracker` object). The values are averaged
along all directions before calculation.

Each correlation function has its own weight as described in
Kirill M. Gerke and Marina V. Karsanina 2015 EPL 111 56002
"""
function euclid_mean_weighted(data1 :: CorrelationTracker,
                              data2 :: CorrelationTracker)
    @assert tracked_data(data1) == tracked_data(data2)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_mean(data1 |> data,
                                       data2 |> data)
                   for data in tracked_data(data1))
    f(:: Any, :: Any) = sumcosts(data1, data2, euclid_mean, weights)
    return f
end

"""
    euclid_directional_weighted(data1 :: CorrelationTracker, data2 :: CorrelationTracker)

Calculate euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `CorrelationTracker` object). The values calculated
along different directions are not averaged and treated separately.

Each correlation function has its own weight as described in
Kirill M. Gerke and Marina V. Karsanina 2015 EPL 111 56002
"""
function euclid_directional_weighted(data1 :: CorrelationTracker,
                                     data2 :: CorrelationTracker)
    @assert tracked_data(data1) == tracked_data(data2)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_directional(data1 |> data,
                                              data2 |> data)
                   for data in tracked_data(data1))
    f(:: Any, :: Any) = sumcosts(data1, data2, euclid_directional, weights)
    return f
end
