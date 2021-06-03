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
    aarts_korst_cooldown(steps = 15, λ = 0.01)

Make the Aarts-Korst cooldown schedule. The temperature is changed
each `steps` steps of annealing algorithm based on standard deviation
of cost function and a parameter `λ`.
"""
function aarts_korst_cooldown(steps :: Integer = 15, λ :: Float64 = 0.01)
    costs = Vector{Float64}(undef, steps)
    counter = 1

    function f(T :: Float64, cost :: Float64)
        costs[counter] = cost
        counter = mod1(counter + 1, steps)

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
    frost_heinemann_cooldown(minimal_steps = 250, λ = 0.01)

Make the Frost-Heinemann cooldown schedule. The temperature is changed
each `steps` steps of annealing algorithm or until some "target"
energy is reached which is based on standard deviation of cost
function and a parameter `λ`.
"""
function frost_heineman_cooldown(minimal_steps :: Integer = 250, λ :: Float64 = 0.01)
    costs = Vector{Float64}(undef, 0)
    sizehint!(costs, minimal_steps)
    target_prev = 0.0
    target = Inf

    function f(T :: Float64, cost :: Float64)
        push!(costs, cost)
        newT = T

        μ = mean(costs)
        if length(costs) == minimal_steps && μ < target
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
