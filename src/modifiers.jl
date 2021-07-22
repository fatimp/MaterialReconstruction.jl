# Modifier types

"""
`AbstractModifier` type is responsible for how configuration update
steps are performed during an annealing process.
"""
abstract type AbstractModifier end

"""
    RandomSwapper()

Create `RandomSwapper` modifier which will swap two random voxels of
different phases during an annealing step.

See also: [`RandomFlipper`](@ref), [`AbstractModifier`](@ref).
"""
struct RandomSwapper <: AbstractModifier end

"""
    RandomFlipper()

Create `RandomFlipper` modifier which will flip phase of a random
voxel in a two-phase system during an annealing step.

See also: [`RandomSwapper`](@ref), [`AbstractModifier`](@ref).
"""
struct RandomFlipper <: AbstractModifier end

"""
    InterfaceSwapper()

Create `InterfaceSwapper` modifier which will swap voxels of different
phases which lie at random points of the interface.

See also: [`RandomSwapper`](@ref), [`RandomFlipper`](@ref),
[`InterfaceFlipper`](@ref), [`AbstractModifier`](@ref).
"""
struct InterfaceSwapper <: AbstractModifier end

"""
    InterfaceFlipper()

Create `InterfaceFlipper` modifier which will flip phase of a random
voxel on an interface of a two-phase system during an annealing step.

See also: [`RandomSwapper`](@ref), [`RandomFlipper`](@ref),
[`InterfaceSwapper`](@ref), [`AbstractModifier`](@ref).
"""
struct InterfaceFlipper <: AbstractModifier end

struct ImprovedRandomFlipper <: AbstractModifier
    histogram :: Vector{Int}
    α         :: Float64
end

# Methods

function modify!(tracker :: AbstractArray, :: RandomSwapper)
    indices = CartesianIndices(tracker)
    index1 = rand(indices)

    while true
        # Try to find element with a different phase
        index2 = rand(indices)
        if tracker[index1] != tracker[index2]
            val1, val2 = tracker[index1], tracker[index2]
            token1 = update_corrfns!(tracker, val2, index1)
            token2 = update_corrfns!(tracker, val1, index2)
            return token1, token2
        end
    end
end

function modify!(tracker :: AbstractArray, :: RandomFlipper)
    indices = CartesianIndices(tracker)
    index = rand(indices)

    return update_corrfns!(tracker, 1 - tracker[index], index)
end

function modify!(tracker :: AbstractArray, :: InterfaceSwapper)
    indices = CartesianIndices(tracker)

    while true
        rndidx = rand(indices)
        ray1 = RandomLineIterator(rndidx)
        ray2 = RandomLineIterator(rndidx)
        ray2 = zip(ray2, drop(ray2, 1))

        start_phase = tracker[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, x) && tracker[x] == start_phase
        index1 = dropwhile(cont, ray1) |> first
        index2 = dropwhile(x -> cont(x[2]), ray2) |> first |> first

        if checkbounds(Bool, tracker, index1) &&
           checkbounds(Bool, tracker, index2)
            val1, val2 = tracker[index1], tracker[index2]
            token1 = update_corrfns!(tracker, val2, index1)
            token2 = update_corrfns!(tracker, val1, index2)
            return token1, token2
        end
    end
end

function modify!(tracker :: AbstractArray, :: InterfaceFlipper)
    indices = CartesianIndices(tracker)

    while true
        rndidx = rand(indices)
        ray = RandomLineIterator(rndidx)
        start_phase = tracker[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, x) && tracker[x] == start_phase
        index = dropwhile(cont, ray) |> first

        if checkbounds(Bool, tracker, index)
            return update_corrfns!(tracker, 1 - tracker[index], index)
        end
    end
end

# Here and later: DPN = different phase neighbors
function count_dpn(array :: AbstractArray{T, N},
                   idx   :: CartesianIndex{N}) where {T, N}
    indices = CartesianIndices(array)
    fidx, lidx = first(indices), last(indices)
    uidx = oneunit(idx)
    elt = array[idx]
    return sum(array[idx] != elt for idx in max(idx - uidx, fidx):min(idx + uidx, lidx))
end

function dpn_histogram(array :: AbstractArray)
    dims = ndims(array)
    hist = zeros(Int, 3^dims)

    for idx in CartesianIndices(array)
        different = count_dpn(array, idx)
        hist[different+1] += 1
    end

    return hist
end

function update_pdf(histogram :: Vector{Int}, α :: Float64)
    prob = [p*α^n for (p, n) in zip(histogram, countfrom(0))]
    prob = prob ./ sum(prob)
    return prob
end

function neighbors_for_update(histogram :: Vector{Int}, α :: Float64)
    prob = update_pdf(histogram, α)
    for idx in 2:length(prob)
        prob[idx] += prob[idx-1]
    end

    r = rand(Float64)
    return findfirst(x -> r < x, prob) - 1
end

@doc raw"""
    ImprovedRandomFlipper(tracker, α)

Create `ImprovedRandomFlipper` modifier which will flip phase of a
random voxel with $n$ different phase neighbors with probability 
$p \propto \alpha^n N(n)$ where $N(n)$ is a number of pixels having
$n$ different phase neighbors in the system.

See also: [`RandomSwapper`](@ref), [`RandomFlipper`](@ref),
[`InterfaceSwapper`](@ref), [`AbstractModifier`](@ref).
"""
ImprovedRandomFlipper(tracker :: CorrelationTracker, α :: Float64) =
    ImprovedRandomFlipper(dpn_histogram(tracker), α)

function flip_and_update_histogram!(tracker  :: AbstractArray{T, N},
                                    modifier :: ImprovedRandomFlipper,
                                    index    :: CartesianIndex{N}) where {T, N}
    histogram = modifier.histogram

    indices = CartesianIndices(tracker)
    fidx, lidx = first(indices), last(indices)
    uidx = oneunit(fidx)

    for idx in max(index - uidx, fidx):min(index + uidx, lidx)
        n_neighbor = count_dpn(tracker, idx)
        histogram[n_neighbor + 1] -= 1
    end

    tracker[index] = 1 - tracker[index]

    for idx in max(index - uidx, fidx):min(index + uidx, lidx)
        n_neighbor = count_dpn(tracker, idx)
        histogram[n_neighbor + 1] += 1
    end

    return index
end

function modify!(tracker :: AbstractArray, modifier :: ImprovedRandomFlipper)
    histogram = modifier.histogram
    α         = modifier.α

    n = neighbors_for_update(histogram, α)
    indices = CartesianIndices(tracker)

    while true
        index = rand(indices)
        if n == count_dpn(tracker, index)
            @assert histogram[n+1] > 0
            return flip_and_update_histogram!(tracker, modifier, index)
        end
    end
end

reject!(tracker :: AbstractArray, modifier :: ImprovedRandomFlipper, state :: CartesianIndex) =
    flip_and_update_histogram!(tracker, modifier, state)

# ? FIXME: Maybe invent some trait here, like ModifierType?
function reject!(tracker :: AbstractArray,
                 _       :: AbstractModifier,
                 state   :: NTuple{2, RollbackToken})
    rollback!(tracker, state[2])
    rollback!(tracker, state[1])
    return nothing
end

reject!(tracker :: AbstractArray,
        _       :: AbstractModifier,
        state   :: RollbackToken) = rollback!(tracker, state)
