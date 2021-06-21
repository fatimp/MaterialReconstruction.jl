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

"""
    C2UpdateProxy(modifier, boxsize = 5)

Create proxy modifier which behaves like `modifier` but also estimates
a number of full c2 updates as a function of number of annealing
steps. Bigger values of `boxsize` result in more accurate estimation
but work slower.
"""
struct C2UpdateProxy{T} <: AbstractModifier
    modifier :: T
    updates  :: Vector{Int}
    boxsize  :: Int
    soft     :: Bool
end

function C2UpdateProxy(modifier :: T;
                       boxsize  :: Int  = 5,
                       soft     :: Bool = false) where T <: AbstractModifier
    updates = zeros(Int, 1)
    sizehint!(updates, 1000000)
    return C2UpdateProxy{T}(modifier, updates, boxsize, soft)
end

# Methods

function modify!(tracker :: CorrelationTracker, :: RandomSwapper)
    indices = CartesianIndices(tracker)
    index1 = rand(indices)

    while true
        # Try to find element with a different phase
        index2 = rand(indices)
        if tracker[index1] != tracker[index2]
            tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
            return index1, index2
        end
    end
end

function modify!(tracker :: CorrelationTracker, :: RandomFlipper)
    indices = CartesianIndices(tracker)
    index = rand(indices)

    tracker[index] = 1 - tracker[index]
    return index
end

function modify!(tracker :: CorrelationTracker, :: InterfaceSwapper)
    indices = CartesianIndices(tracker)

    while true
        rndidx = rand(indices)
        ray1 = RandomLineIterator(rndidx)
        ray2 = RandomLineIterator(rndidx)
        ray2 = zip(ray2, drop(ray2, 1))

        start_phase = tracker[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, x) && tracker[x] == start_phase
        idx1 = dropwhile(cont, ray1) |> first
        idx2 = dropwhile(x -> cont(x[2]), ray2) |> first |> first

        if checkbounds(Bool, tracker, idx1) && checkbounds(Bool, tracker, idx2)
            tracker[idx1], tracker[idx2] = tracker[idx2], tracker[idx1]
            return idx1, idx2
        end
    end
end

function modify!(tracker :: CorrelationTracker, :: InterfaceFlipper)
    indices = CartesianIndices(tracker)

    while true
        rndidx = rand(indices)
        ray = RandomLineIterator(rndidx)
        start_phase = tracker[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, x) && tracker[x] == start_phase
        index = dropwhile(cont, ray) |> first

        if checkbounds(Bool, tracker, index)
            tracker[index] = 1 - tracker[index]
            return index
        end
    end
end

function c2_counter_increment(tracker :: CorrelationTracker,
                              idx     :: CartesianIndex{N},
                              boxsize :: Int,
                              soft    :: Bool) where N
    # Take a box with a size of 2boxsize+1 and idx in center
    box = CartesianIndices(Tuple((x - boxsize):(x + boxsize) for x in Tuple(idx)))
    slice = CircularArray(tracker)[box] |> Array

    # Segmentation after update
    segments1 = label_components(slice)

    # Roll back and get segmentation before update
    idx = CartesianIndex(Tuple(fill(boxsize + 1, N)))
    slice[idx] = 1 - slice[idx]
    segments2 = label_components(slice)

    # Compare segmentation before and after
    nonequal = segments1 .!= segments2

    if soft
        # Soft updates: Check if changes in segmentation go beyond the
        # box.
        box_inner = CartesianIndices(Tuple(fill(2:2boxsize, N)))
        nonequal_inner = nonequal[box_inner]
        return (sum(nonequal) == sum(nonequal_inner)) ? 0 : 1
    else
        # Hard updates: Check if segmentation stays the same
        return (sum(nonequal) == 1) ? 0 : 1
    end
end

function modify!(tracker :: CorrelationTracker, modifier :: C2UpdateProxy)
    boxsize = modifier.boxsize
    updates = modifier.updates
    soft    = modifier.soft

    idx :: CartesianIndex = modify!(tracker, modifier.modifier)
    push!(updates, updates[end] + c2_counter_increment(tracker, idx, boxsize, soft))
    return idx
end

# ? FIXME: Maybe invent some trait here, like ModifierType?
function rollback!(tracker :: CorrelationTracker,
                           :: AbstractModifier,
                   state   :: Tuple{CartesianIndex, CartesianIndex})
    index1, index2 = state
    tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
end

function rollback!(tracker :: CorrelationTracker,
                           :: AbstractModifier,
                   state   :: CartesianIndex)
    tracker[state] = 1 - tracker[state]
end
