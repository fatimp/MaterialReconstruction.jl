# Different cooldown schedules. Internally, a cooldown function
# accepts the current temperature as its first argument and the
# current energy as its second argument and return the new temperature
# of the system.

@doc raw"""
    exponential_cooldown(λ :: Float64)

Create an exponential cooldown schedule. Temperature on step $n + 1$
is calculated as $T_{n + 1} = \lambda T_{n}$.
"""
exponential_cooldown(λ :: Float64 = 0.999999) = (T :: Float64, :: Float64) -> λ * T

"""
    aarts_korst_cooldown(;n = 15, λ = 0.01)

Make the Aarts-Korst cooldown schedule. The temperature decreases
each `n` steps of annealing algorithm. The higer `λ` is the faster
decreases the temperature.

For more information, see Aarts, E.H.L. and Korst, J.H.M. (1989)
Simulated Annealing and Boltzmann Machines: A Stochastic Approach to
Combinatorial Optimization and Neural Computing. John Wiley & Sons,
Chichester.
"""
function aarts_korst_cooldown(;n :: Integer = 15, λ :: Float64 = 0.01)
    costs = Vector{Float64}(undef, n)
    counter = 1

    function f(T :: Float64, cost :: Float64)
        costs[counter] = cost
        counter = mod1(counter + 1, n)

        newT = T
        if counter == 1
            σ = std(costs)
            newT = T * σ / (σ + λ*T)
        end

        return newT
    end

    return f
end

"""
    frost_heineman_cooldown(;n = 250, λ = 0.01)

Make the Frost-Heineman cooldown schedule. The temperature decreases
each `n` steps of annealing algorithm or until some "target" energy is
reached. The target energy is based on standard deviation of the cost
function and a parameter `λ`.

For more information, see R. Frost, P. Heineman "Simulated Annealing:
A Heuristic for Parallel Stochastic Optimization" (1997)
"""
function frost_heineman_cooldown(;n :: Integer = 250, λ :: Float64 = 0.01)
    costs = Vector{Float64}(undef, 0)
    sizehint!(costs, n)
    target = Inf

    function f(T :: Float64, cost :: Float64)
        push!(costs, cost)
        newT = T

        μ = mean(costs)
        if length(costs) == n && μ < target
            σ = std(costs)
            empty!(costs)

            target_prev = isinf(target) ? μ : target
            target = μ - λ*σ
            newT = T + (target - target_prev) * (T/σ)^2
        end

        return newT
    end

    return f
end
