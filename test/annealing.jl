function test_annealing(target, init, cost, modifier, cooldown)
    system = CorrelationTracker(target |> init, target)
    furnace = Furnace(system, target; T0 = 0.0)
    c0 = euclid_directional(system, target)

    for i in 1:10000
        furnace = annealing_step(furnace;
                                 cost     = cost,
                                 modifier = modifier,
                                 cooldown = cooldown)
    end

    c1 = euclid_directional(system, target)

    @test c1 < c0
end

@testset "300x300 image" begin
    array = Utilities.read_cuboid("image3d.json")[:,:,1]
    target = CorrelationTracker(array;
                                directions = [:x, :y, :xy, :yx],
                                periodic   = true)

    # Need to find starting costs for this
    #modifiers = (InterfaceFlipper(), InterfaceSwapper(),
    #             RandomSwapper(),    RandomFlipper())
    #costs = (euclid_mean, euclid_directional)
    #initializers = (initialize_random, initialize_spheres)
    #cooldowns = (exponential_cooldown(), aarts_korst_cooldown(), frost_heineman_cooldown())

    #alltogether = product(modifiers, costs, initializers, cooldowns)

    alltogether = product((Flipper(InterfaceSampler()), Swapper(InterfaceSampler())),
                          (euclid_directional, euclid_mean),
                          (initialize_spheres,),
                          (aarts_korst_cooldown(), exponential_cooldown()))

    foreach(alltogether) do prod
        modifier, cost, initializer, cooldown = prod
        test_annealing(target, initializer, cost, modifier, cooldown)
    end
end
