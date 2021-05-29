struct Furnace{T, N}
    system      :: CorrelationTracker{T, N}
    target      :: CorrelationTracker{T, N}
    temperature :: Float64
    steps       :: Int
    rejected    :: Int
    # This is somewhat misleading. It counts accepted permutations
    # which result in an increase of metric.
    accepted    :: Int
end

# Initialize a system for annealing with random data, but preserving
# porosity of the target
function initialize_system(target :: CorrelationTracker{T, N},
                           shape  :: NTuple{N, Int}) where {T, N}
    # The porosity may be already calculated as $S^1_2(0)$ or
    # $L^1_2(0)$, but what if we do not use these functions for
    # reconstruction? Better to recalculate it again, it's fast.
    porosity = count(x -> x == 1, target) / length(target)

    # Only two-phase systems are supported by now. So to begin the
    # reconstruction we start with a system with all elements
    # belonging to the phase "zero".
    array = zeros(T, shape)

    # Determine the number of elements with phase "one" to keep the
    # porosity of the original
    n = porosity * prod(shape)

    rndidx(shape) = CartesianIndex((rand(1:s) for s in shape)...)

    # Populate the array
    while n > 0
        idx = rndidx(shape)
        # If we have a free place, put phase "one" to it
        if array[idx] == 0
            array[idx] = 1
            n -= 1
        end
    end

    # KLUDGE: Maybe add a function to return directions and length?
    dirs = []
    len  = 0
    for data in tracked_data(target)
        dirs = target |> data |> directions
        len  = target |> data |> length
        break
    end
    
    # Create CorrelationTracker for the system being annealed
    system = CorrelationTracker{T, N}(array;
                                      tracking   = target |> tracked_data |> collect,
                                      periodic   = target.periodic,
                                      directions = dirs,
                                      len        = len)
end

Furnace(target :: CorrelationTracker{T, N},
        shape  :: NTuple{N, Int};
        T0     :: Float64) where {T, N} =
            Furnace(initialize_system(target, shape), target, T0, 0, 0, 0)

Base.show(io :: IO, furnace :: Furnace) = begin
    print(io, "Furnace with system of dimensions $(size(furnace.system)), T=$(furnace.temperature) ")
    print(io, "steps=$(furnace.steps), rejected=$(furnace.rejected), accepted=$(furnace.accepted)")
end
