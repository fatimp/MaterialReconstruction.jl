using CorrelationTrackers
using CorrelationFunctions
using MaterialReconstruction
using XUnit
using Base.Iterators: product

@testset "Test simulated annealing" begin include("annealing.jl") end
