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

Create a sampler which picks samples randomly from an interface
between phases.

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
    DPNSampler(array, α)

Create a sampler which picks a random sample with $n$ different phase
neighbors with probability  $p \propto \alpha^n N(n)$ where $N(n)$ is
a number of pixels having $n$ different phase neighbors in the system.

See also: [`UniformSampler`](@ref), [`InterfaceSampler`](@ref),
[`AbstractSampler`](@ref).
"""
DPNSampler(array :: AbstractArray, α :: Float64) =
    DPNSampler(dpn_histogram(array), α)

"""
`AbstractModifier` type is responsible for how configuration update
steps are performed during an annealing process.

See also: [`Flipper`](@ref), [`Swapper`](@ref).
"""
abstract type AbstractModifier end

"""
    Flipper(sampler)

Create a modifier which changes phase of a sample taken from a binary
array.

See also: [`AbstractModifier`](@ref), [`Swapper`](@ref).
"""
struct Flipper{T <: AbstractSampler} <: AbstractModifier
    sampler :: T
end

"""
    Swapper(sampler)

Create a modifier which takes two samples having different phases and
swaps phases of these two samples.

See also: [`Flipper`](@ref), [`AbstractModifier`](@ref).
"""
struct Swapper{T <: AbstractSampler} <: AbstractModifier
    sampler :: T
end

# Interface

"""
    sample(array, sampler)

Get a random sample from the array. This function must return an
index, e.g. `CartesianIndex`.
"""
function sample end

"""
    update_pre!(array, index, sampler)

Update the state of sampler before the value `array[index]` is
changed. Should only be implemented for stateful samplers.
"""
function update_pre! end

"""
    update_post!(array, index, sampler)

Update the state of sampler after the value `array[index]` is
changed. Should only be implemented for stateful samplers.
"""
function update_post! end

"""
    modify!(array, modifier)

Randomly modify the array. This function must return a state used in
`reject!` to bring the array to a previous state if this modification
is rejected.
"""
function modify! end

"""
    reject!(array, modifier, state)

Bring the array back to the previous state.
"""
function reject! end

# Methods

# Stateless samples do not have to implement this
update_pre!( :: Any, :: Any, :: AbstractSampler) = nothing
update_post!(:: Any, :: Any, :: AbstractSampler) = nothing

sample(array :: AbstractArray, :: UniformSampler) =
    array |> CartesianIndices |> rand

function sample(array :: AbstractArray, :: InterfaceSampler)
    indices = CartesianIndices(array)

    while true
        rndidx = rand(indices)
        ray = RandomLineIterator(rndidx)
        start_phase = array[rndidx]

        cont(x :: CartesianIndex) =
            checkbounds(Bool, array, x) && array[x] == start_phase
        index = dropwhile(cont, ray) |> first

        if checkbounds(Bool, array, index)
            return index
        end
    end
end


# That horrible stateful sampler
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

function update_histogram!(array   :: AbstractArray{T, N},
                           index   :: CartesianIndex{N},
                           sampler :: DPNSampler,
                           fn      :: Function) where {T, N}
    histogram = sampler.histogram

    indices = CartesianIndices(array)
    fidx, lidx = first(indices), last(indices)
    uidx = oneunit(fidx)

    for idx in max(index - uidx, fidx):min(index + uidx, lidx)
        n_neighbor = count_dpn(array, idx)
        histogram[n_neighbor + 1] = fn(histogram[n_neighbor + 1])
    end
end

update_pre!(array   :: AbstractArray{T, N},
            index   :: CartesianIndex{N},
            sampler :: DPNSampler) where {T, N} =
                update_histogram!(array, index, sampler, x -> x - 1)

update_post!(array   :: AbstractArray{T, N},
             index   :: CartesianIndex{N},
             sampler :: DPNSampler) where {T, N} =
                 update_histogram!(array, index, sampler, x -> x + 1)

function sample(array :: AbstractArray, sampler :: DPNSampler)
    histogram = sampler.histogram
    α         = sampler.α

    n = neighbors_for_update(histogram, α)
    indices = CartesianIndices(array)
    @assert histogram[n + 1] > 0

    while true
        index = rand(indices)
        if n == count_dpn(array, index)
            return index
        end
    end
end

# Update correlation functions and sampler state
function update_all!(array   :: AbstractArray{T, N},
                     sampler :: AbstractSampler,
                     index   :: CartesianIndex{N},
                     val) where {T, N}
    update_pre!(array, index, sampler)
    token = update_corrfns!(array, val, index)
    update_post!(array, index, sampler)

    return token
end

# Rollback correlation functions and sampler state
function reject_all!(array   :: AbstractArray{T, N},
                     sampler :: AbstractSampler,
                     token   :: AbstractRollbackToken) where {T, N}
    update_pre!(array, token.index, sampler)
    rollback!(array, token)
    update_post!(array, token.index, sampler)

    return nothing
end

function modify!(array :: AbstractArray, modifier :: Flipper)
    sampler = modifier.sampler
    index   = sample(array, modifier.sampler)

    return update_all!(array, sampler, index, 1 - array[index])
end

reject!(array    :: AbstractArray,
        modifier :: Flipper,
        state    :: AbstractRollbackToken) =
            reject_all!(array, modifier.sampler, state)

function modify!(array :: AbstractArray, modifier :: Swapper)
    sampler = modifier.sampler
    index1 = sample(array, sampler)
    val1 = array[index1]

    while true
        index2 = sample(array, modifier.sampler)
        val2 = array[index2]
        if val1 ≠ val2
            token1 = update_all!(array, sampler, index1, val2)
            token2 = update_all!(array, sampler, index2, val1)
            return token1, token2
        end
    end
end

function reject!(array    :: AbstractArray,
                 modifier :: Swapper,
                 state    :: NTuple{2, AbstractRollbackToken})
    token1, token2 = state
    sampler = modifier.sampler

    reject_all!(array, sampler, token2)
    reject_all!(array, sampler, token1)
end
