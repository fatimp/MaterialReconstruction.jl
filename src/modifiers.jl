# Modifier and sampler types

"""
`AbstractSampler` type is responsible for picking samples (voxels or
pixels) during annealing process.

See also: [`UniformSampler`](@ref), [`InterfaceSampler`](@ref),
[`DPNSampler`](@ref).
"""
abstract type AbstractSampler end

"""
    UniformSampler()

Create a sampler which picks samples randomly with uniform
distribution.

See also: [`AbstractSampler`](@ref), [`InterfaceSampler`](@ref),
[`DPNSampler`](@ref).
"""
struct UniformSampler <: AbstractSampler end

"""
    InterfaceSampler()

Create a sampler which picks samples randomly from two-phase
interface.

See also: [`UniformSampler`](@ref), [`AbstractSampler`](@ref),
[`DPNSampler`](@ref).
"""
struct InterfaceSampler <: AbstractSampler end

# Here and later: DPN = different phase neighbors
struct DPNSampler <: AbstractSampler
    histogram :: Vector{Int}
    α         :: Float64
end

@doc raw"""
    DPNSampler(tracker, α)

Create a sampler which picks a random sampler with $n$ different phase
neighbors with probability  $p \propto \alpha^n N(n)$ where $N(n)$ is
a number of pixels having $n$ different phase neighbors in the
system.

See also: [`UniformSampler`](@ref), [`InterfaceSampler`](@ref),
[`AbstractSampler`](@ref).
"""
DPNSampler(tracker :: CorrelationTracker, α :: Float64) =
    DPNSampler(dpn_histogram(tracker), α)

"""
`AbstractModifier` type is responsible for how configuration update
steps are performed during an annealing process.
"""
abstract type AbstractModifier end

"""
    Flipper(sampler)

Create a modifier which flips phase of a sample during annealing.
"""
struct Flipper{T <: AbstractSampler} <: AbstractModifier
    sampler :: T
end

"""
    Swapper(sampler)

Create a modifier which swaps phases of two samples during annealing.
"""
struct Swapper{T <: AbstractSampler} <: AbstractModifier
    sampler :: T
end

# Methods

sample(tracker :: AbstractArray, :: UniformSampler) =
    tracker |> CartesianIndices |> rand

function sample(tracker :: AbstractArray, :: InterfaceSampler)
    indices = CartesianIndices(tracker)

    while true
        rndidx = rand(indices)
        ray = RandomLineIterator(rndidx)
        start_phase = tracker[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, tracker, x) && tracker[x] == start_phase
        index = dropwhile(cont, ray) |> first

        if checkbounds(Bool, tracker, index)
            return index
        end
    end
end

function modify!(tracker :: AbstractArray, modifier :: Flipper)
    index = sample(tracker, modifier.sampler)
    return update_corrfns!(tracker, 1 - tracker[index], index)
end

reject!(tracker :: AbstractArray,
        _       :: Flipper,
        state   :: RollbackToken) = rollback!(tracker, state)

function modify!(tracker :: AbstractArray, modifier :: Swapper)
    index1 = sample(tracker, modifier.sampler)
    val1 = tracker[index1]

    while true
        index2 = sample(tracker, modifier.sampler)
        val2 = tracker[index2]
        if val1 ≠ val2
            token1 = update_corrfns!(tracker, val2, index1)
            token2 = update_corrfns!(tracker, val1, index2)
            return token1, token2
        end
    end
end

function reject!(tracker :: AbstractArray,
                 _       :: Swapper,
                 state   :: NTuple{2, RollbackToken})
    rollback!(tracker, state[2])
    rollback!(tracker, state[1])
    return nothing
end

# function count_dpn(array :: AbstractArray{T, N},
#                    idx   :: CartesianIndex{N}) where {T, N}
#     indices = CartesianIndices(array)
#     fidx, lidx = first(indices), last(indices)
#     uidx = oneunit(idx)
#     elt = array[idx]
#     return sum(array[idx] != elt for idx in max(idx - uidx, fidx):min(idx + uidx, lidx))
# end

# function dpn_histogram(array :: AbstractArray)
#     dims = ndims(array)
#     hist = zeros(Int, 3^dims)
#
#     for idx in CartesianIndices(array)
#         different = count_dpn(array, idx)
#         hist[different+1] += 1
#     end
#
#     return hist
# end

# function update_pdf(histogram :: Vector{Int}, α :: Float64)
#     prob = [p*α^n for (p, n) in zip(histogram, countfrom(0))]
#     prob = prob ./ sum(prob)
#     return prob
# end

# function neighbors_for_update(histogram :: Vector{Int}, α :: Float64)
#     prob = update_pdf(histogram, α)
#     for idx in 2:length(prob)
#         prob[idx] += prob[idx-1]
#     end
#
#     r = rand(Float64)
#     return findfirst(x -> r < x, prob) - 1
# end

# function flip_and_update_histogram!(tracker  :: AbstractArray{T, N},
#                                     modifier :: ImprovedRandomFlipper,
#                                     index    :: CartesianIndex{N}) where {T, N}
#     histogram = modifier.histogram
#
#     indices = CartesianIndices(tracker)
#     fidx, lidx = first(indices), last(indices)
#     uidx = oneunit(fidx)
#
#     for idx in max(index - uidx, fidx):min(index + uidx, lidx)
#         n_neighbor = count_dpn(tracker, idx)
#         histogram[n_neighbor + 1] -= 1
#     end
#
#     tracker[index] = 1 - tracker[index]
#
#     for idx in max(index - uidx, fidx):min(index + uidx, lidx)
#         n_neighbor = count_dpn(tracker, idx)
#         histogram[n_neighbor + 1] += 1
#     end
#
#     return index
# end

# function modify!(tracker :: AbstractArray, modifier :: ImprovedRandomFlipper)
#     histogram = modifier.histogram
#     α         = modifier.α
#
#     n = neighbors_for_update(histogram, α)
#     indices = CartesianIndices(tracker)
#
#     while true
#         index = rand(indices)
#         if n == count_dpn(tracker, index)
#             @assert histogram[n+1] > 0
#             return flip_and_update_histogram!(tracker, modifier, index)
#         end
#     end
# end

# reject!(tracker :: AbstractArray, modifier :: ImprovedRandomFlipper, state :: CartesianIndex) =
#     flip_and_update_histogram!(tracker, modifier, state)
