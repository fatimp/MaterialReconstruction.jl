# Line "drawing" iterator
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
