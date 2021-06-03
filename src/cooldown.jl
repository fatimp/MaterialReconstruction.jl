# Different cooldown schedules. Internally, a cooldown function
# accepts the current temperature as its first argument and the
# current energy as its second argument and return the new temperature
# of the system.

@doc raw"""
    exponential_cooldown(λ :: Float64)

Make an exponential cooldown function. Temperature on step $n + 1$ is
calculated as $T_{n + 1} = \lambda T_{n}$.
"""
exponential_cooldown(λ :: Float64 = 0.999999) = (T :: Float64, :: Float64) -> λ * T
