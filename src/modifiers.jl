"""
`AbstractModifier` type is responsible for how configuration update
steps are performed during an annealing process.
"""
abstract type AbstractModifier end

# ? FIXME: Maybe invent some trait here, like ModifierType?
function rollback!(modifier :: AbstractModifier,
                   state    :: Tuple{CartesianIndex, CartesianIndex})
    tracker = modifier.tracker
    index1, index2 = state
    tracker[index1], tracker[index2] = tracker[index2], tracker[index1]
end

function rollback!(modifier :: AbstractModifier,
                   state    :: CartesianIndex)
    tracker = modifier.tracker
    tracker[state] = 1 - tracker[state]
end

"""
    RandomSwapper(tracker :: CorrelationTracker)

Create `RandomSwapper` modifier which will swap two random voxels of
different phases during an annealing step.

See also: [`RandomFlipper`](@ref), [`AbstractModifier`](@ref).
"""
struct RandomSwapper <: AbstractModifier
    tracker :: CorrelationTracker
end

function modify!(modifier :: RandomSwapper)
    tracker = modifier.tracker
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

"""
    RandomFlipper(tracker :: CorrelationTracker)

Create `RandomFlipper` modifier which will flip phase of a random
voxel in a two-phase system during an annealing step.

See also: [`RandomSwapper`](@ref), [`AbstractModifier`](@ref).
"""
struct RandomFlipper <: AbstractModifier
    tracker :: CorrelationTracker
end

function modify!(modifier :: RandomFlipper)
    tracker = modifier.tracker
    shape = size(tracker)
    index = random_index(shape)

    tracker[index] = 1 - tracker[index]
    return index
end

# Line "drawing iterator"
struct LineIterator{N}
    start :: NTuple{N, Int}
    ϕ     :: Float64
    θ     :: Float64
end

RandomLineIterator(start) =
    LineIterator(start, 2π*rand(Float64), π*rand(Float64) - π/2)

Base.IteratorSize(::LineIterator) = Base.IsInfinite()
Base.iterate(iter :: LineIterator) = CartesianIndex(iter.start), 0.1

function Base.iterate(iter :: LineIterator{2}, r :: Float64)
    r    = r + sqrt(2)
    x, y = iter.start
    ϕ    = iter.ϕ

    xn = x + r*cos(ϕ) |> floor |> Int
    yn = y + r*sin(ϕ) |> floor |> Int
    return CartesianIndex(xn, yn), r
end

function Base.iterate(iter :: LineIterator{3}, r :: Float64)
    r       = r + sqrt(3)
    x, y, z = iter.start
    ϕ       = iter.ϕ
    θ       = iter.θ

    xn = x + r*cos(θ)*cos(ϕ) |> floor |> Int
    yn = y + r*cos(θ)*sin(ϕ) |> floor |> Int
    zn = z + r*sin(θ)        |> floor |> Int
    return CartesianIndex(xn, yn, zn), r
end

"""
    InterfaceSwapper(tracker :: CorrelationTracker)

Create `InterfaceSwapper` modifier which will swap voxels of different
phases which lie at random points of the interface.

See also: [`RandomSwapper`](@ref), [`RandomFlipper`](@ref),
[`AbstractModifier`](@ref).
"""
struct InterfaceSwapper <: AbstractModifier
    tracker :: CorrelationTracker
end

function modify!(modifier :: InterfaceSwapper)
    tracker = modifier.tracker
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
