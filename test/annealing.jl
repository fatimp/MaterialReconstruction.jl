function test_annealing(target, init, cost, modifier, cooldown)
    system = init(target, (1000, 1000))
    furnace = Furnace(system, target; T0 = 7e-5)
    c = cost(system, target)
    c0 = euclid_directional(system, target)

    for i in 1:5000
        furnace = annealing_step(furnace;
                                 cost     = c,
                                 modifier = modifier,
                                 cooldown = cooldown)
    end

    c1 = euclid_directional(system, target)

    @test c1 < c0
end

@testset "Value noise original" begin
    array = Int8.([value_noise(x/50, y/100, 0.0, 4, 43565) for x in 1:1000, y in 1:1000] .< 0.5)
    target = CorrelationTracker{Int8, 2}(array; directions = [:x, :y, :xy_main, :xy_anti])

    test_annealing(target, initialize_spheres, čapek_cost, InterfaceFlipper(), aarts_korst_cooldown())
end
