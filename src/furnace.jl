struct Furnace{A1 <: AbstractArray, A2 <: AbstractArray}
    system      :: A1
    target      :: A2
    temperature :: Float64
    steps       :: Int
    rejected    :: Int
    # This is somewhat misleading. It counts accepted permutations
    # which result in an increase of metric.
    accepted    :: Int
end

"""
    Furnace(system :: AbstractArray, target :: AbstractArray; T0)

Initialize a furnace (an object which is used in annealing
process). `system` is an array being reconstructed to be similar to
`target`. `T0` is an initial temperature of a furnace.
"""
Furnace(system :: A1,
        target :: A2;
        T0     :: Float64) where {A1 <: AbstractArray, A2 <: AbstractArray} =
            Furnace(system, target, T0, 0, 0, 0)

Base.show(io :: IO, furnace :: Furnace) = begin
    print(io, "Furnace with system of dimensions $(size(furnace.system)), T=$(furnace.temperature) ")
    print(io, "steps=$(furnace.steps), rejected=$(furnace.rejected), accepted=$(furnace.accepted)")
end
