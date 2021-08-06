metric(x :: AbstractArray, y :: AbstractArray) =
    mapreduce((x, y) -> (x - y)^2, +, x, y)

function sumcosts(tracker1 :: AbstractArray{T, N},
                  tracker2 :: AbstractArray{T, N},
                  costfn   :: Function) where {T, N}
    @assert tracked_data(tracker1) == tracked_data(tracker2)
    return mapreduce(data -> costfn(tracker1 |> data,
                                    tracker2 |> data),
                     +, tracked_data(tracker1))
end

function sumcosts(tracker1 :: AbstractArray{T, N},
                  tracker2 :: AbstractArray{T, N},
                  costfn   :: Function,
                  weights  :: Dict{AbstractTracker{T}, Float64}) where {T, N}
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
    euclid_mean(data1 :: AbstractArray, data2 :: AbstractArray)

Calculate squared euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `AbstractArray` object). The values are averaged
along all directions before calculation.
"""
function euclid_mean end

euclid_mean(data1 :: CorrelationData{T},
            data2 :: CorrelationData{T}) where T =
    metric(data1 |> mean, data2 |> mean)

euclid_mean(data1 :: AbstractArray{T, N},
            data2 :: AbstractArray{T, N}) where {T, N} =
                sumcosts(data1, data2, euclid_mean)

"""
    euclid_directional(data1 :: CorrelationData, data2 :: CorrelationData)
    euclid_directional(data1 :: AbstractArray, data2 :: AbstractArray)

Calculate squared euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `AbstractArray` object). The values calculated
along different directions are not averaged and treated separately.
"""
function euclid_directional end

function euclid_directional(data1 :: CorrelationData{T},
                            data2 :: CorrelationData{T}) where T
    @assert directions(data1) == directions(data2)
    sum(metric(data1[dir], data2[dir]) for dir in directions(data1))
end

function euclid_directional(data1 :: AbstractArray{T, N},
                            data2 :: AbstractArray{T, N}) where {T, N}
    sumcosts(data1, data2, euclid_directional)
end

"""
    euclid_mean_weighted(data1 :: AbstractArray, data2 :: AbstractArray)

Calculate squared euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `AbstractArray` object). The values are averaged
along all directions before calculation.

Each correlation function has its own weight as described in
Kirill M. Gerke and Marina V. Karsanina 2015 EPL 111 56002
"""
function euclid_mean_weighted(data1 :: AbstractArray{T, N},
                              data2 :: AbstractArray{T, N}) where {T, N}
    @assert tracked_data(data1) == tracked_data(data2)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_mean(data1 |> data,
                                       data2 |> data)
                   for data in tracked_data(data1))
    f(:: Any, :: Any) = sumcosts(data1, data2, euclid_mean, weights)
    return f
end

"""
    euclid_directional_weighted(data1 :: AbstractArray, data2 :: AbstractArray)

Calculate squared euclidean distance between values of two correlation
functions (represented as `CorrelationData` object) or two systems
(represented as `AbstractArray` object). The values calculated
along different directions are not averaged and treated separately.

Each correlation function has its own weight as described in
Kirill M. Gerke and Marina V. Karsanina 2015 EPL 111 56002
"""
function euclid_directional_weighted(data1 :: AbstractArray{T, N},
                                     data2 :: AbstractArray{T, N}) where {T, N}
    @assert tracked_data(data1) == tracked_data(data2)
    # At first, calculate weights for all correlation functions
    weights = Dict(data => euclid_directional(data1 |> data,
                                              data2 |> data)
                   for data in tracked_data(data1))
    f(:: Any, :: Any) = sumcosts(data1, data2, euclid_directional, weights)
    return f
end

"""
    generalized_čapek_cost(tracker1, tracker2, dict)

Returns a function which calculates the cost based on two-point and
lineal-path functions for solid phase and any other functions present
in `keys(dict)`. `dict` is a dictionary which includes key-value pairs
`AbstractTracker => Float64`. Values of `dict` must be in range
`[0, 1]`. They control an initial contribution of the corresponding
correlation function. The smaller these values are the smaller is
the initial contribution. Contributions of all correlation functions
equivalizes when the cost function becomes small.
"""
function generalized_čapek_cost(data1 :: AbstractArray{T, N},
                                data2 :: AbstractArray{T, N},
                                dict  :: Dict{<:AbstractTracker{T}, Float64}) where {T, N}
    # Initial S2
    s2cost_init  = euclid_directional(Directional.s2(data1, 0),
                                      Directional.s2(data2, 0))
    # Initial L2 for solid phase
    l2scost_init = euclid_directional(Directional.l2(data1, 1),
                                      Directional.l2(data2, 1))

    # Rescale coefficients as in čapek_cost
    sdict = Dict(k => v*(s2cost_init + l2scost_init) for (k, v) in dict)

    local function f(data1 :: AbstractArray, data2 :: AbstractArray)
        # S2
        s2cost  = euclid_directional(Directional.s2(data1, 0),
                                     Directional.s2(data2, 0))
        # L2 for solid phase
        l2scost = euclid_directional(Directional.l2(data1, 1),
                                     Directional.l2(data2, 1))
        # All other cost functions
        other = sum(η/(η + s2cost + l2scost) *
                    euclid_directional(corrfn(data1),
                                       corrfn(data2))
                    for (corrfn, η) in sdict)

        return s2cost + l2scost + other
    end

    return f
end

"""
    čapek_cost(data1 :: AbstractArray, data2 :: AbstractArray, η = 0.6)

Returns a function which calculates the cost based on two-point
correlation function for void phase  and lineal-path function for
solid and void phases where contribution of the latter increases with
time.
"""
čapek_cost(data1 :: AbstractArray{T, N},
           data2 :: AbstractArray{T, N},
           η     :: Float64 = 0.6) where {T, N} =
               generalized_čapek_cost(data1, data2, Dict(T |> zero |> L2Tracker => η))
