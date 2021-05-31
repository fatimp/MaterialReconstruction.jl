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

"""
    Furnace(system :: CorrelationTracker, target :: CorrelationTracker; T0)

Initialize a furnace (an object which is used in annealing
process). `system` is a system being reconstructed. `target` holds
desired values of correlation functions. `T0` is an initial
temperature of a furnace.
"""
Furnace(system :: CorrelationTracker{T, N},
        target :: CorrelationTracker{T, N};
        T0     :: Float64) where {T, N} =
            Furnace(system, target, T0, 0, 0, 0)

Base.show(io :: IO, furnace :: Furnace) = begin
    print(io, "Furnace with system of dimensions $(size(furnace.system)), T=$(furnace.temperature) ")
    print(io, "steps=$(furnace.steps), rejected=$(furnace.rejected), accepted=$(furnace.accepted)")
end
