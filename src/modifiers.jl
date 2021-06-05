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

# Methods

function modify!(tracker :: CorrelationTracker, :: RandomSwapper)
    shape = size(tracker)
    index1 = random_index(shape)

    while true
        # Try to find element with a different phase
        index2 = random_index(shape)
        if tracker[index1] != tracker[index2]
            tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
            return index1, index2
        end
    end
end

function modify!(tracker :: CorrelationTracker, :: RandomFlipper)
    shape = size(tracker)
    index = random_index(shape)

    tracker[index] = 1 - tracker[index]
    return index
end

function modify!(tracker :: CorrelationTracker, :: InterfaceSwapper)
    shape = size(tracker)

    while true
        rndidx = Tuple(rand(1:x) for x in shape)
        ray1 = RandomLineIterator(rndidx)
        ray2 = RandomLineIterator(rndidx)
        ray2 = zip(ray2, drop(ray2, 1))

        start_phase = tracker[rndidx...]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, Tuple(x)...) && tracker[x] == start_phase
        idx1 = dropwhile(cont, ray1) |> first
        idx2 = dropwhile(x -> cont(x[2]), ray2) |> first |> first

        if checkbounds(Bool, tracker, Tuple(idx1)...) &&
           checkbounds(Bool, tracker, Tuple(idx2)...)
            tracker[idx1], tracker[idx2] = tracker[idx2], tracker[idx1]
            return idx1, idx2
        end
    end
end

function modify!(tracker :: CorrelationTracker, :: InterfaceFlipper)
    shape = size(tracker)

    while true
        rndidx = Tuple(rand(1:x) for x in shape)
        ray = RandomLineIterator(rndidx)
        start_phase = tracker[rndidx...]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, Tuple(x)...) && tracker[x] == start_phase
        index = dropwhile(cont, ray) |> first

        if checkbounds(Bool, tracker, Tuple(index)...)
            tracker[index] = 1 - tracker[index]
            return index
        end
    end
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
