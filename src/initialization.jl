"""
    initialize_random(target :: AbstractArray[, shape])

Initialize a system for the annealing procedure. The system is to be
filled with random data in such a way that the porosity of the
`target` is preserved.
"""
function initialize_random(target :: AbstractArray{T, N},
                           shape  :: NTuple{N, Int} = size(target)) where {T, N}
    # The porosity may be already calculated as $S^1_2(0)$ or
    # $L^1_2(0)$, but what if we do not use these functions for
    # reconstruction? Better to recalculate it again, it's fast.
    porosity = count(x -> x == 1, target) / length(target)

    # Only two-phase systems are supported by now. So to begin the
    # reconstruction we start with a system with all elements
    # belonging to the phase "zero".
    array = zeros(T, shape)

    # All indices in array
    indices = CartesianIndices(array)

    # Determine the number of elements with phase "one" to keep the
    # porosity of the original
    n = porosity * prod(shape)

    # Populate the array
    while n > 0
        idx = rand(indices)
        # If we have a free place, put phase "one" to it
        if array[idx] == 0
            array[idx] = 1
            n -= 1
        end
    end

    return array
end

"""
    initialize_spheres(target :: AbstractArray[, shape])

Initialize a system for the annealing procedure. The system is to be
filled with random spheres so its two-point correlation function for
the void phase is as close as possible to the `target`'s one.
"""
function initialize_spheres(target :: AbstractArray{T, N};
                            shape  :: NTuple{N, Int}  = size(target),
                            R0     :: Float64         = 10.0,
                            λ0     :: Float64         = 1e-3) where {T, N}
    R, λ = spheres_parameters(target, R0, λ0)
    return generate_spheres(shape, T, R, λ)
end
