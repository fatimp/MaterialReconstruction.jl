# MaterialReconstruction.jl

`MaterialReconstruction.jl` package can be used to recreate binary arrays from
a set of descriptors known as correlation functions using a technique called
simulated annealing. It works in conjunction with packages
`CorrelationFunctions.jl` and `CorrelationsTrackers.jl`.

This package is quite configurable which means you can decide which correlation
functions to take into account, which cooldown schedule to use etc.

Look at this minimal example to see how `MaterialReconstruction.jl` can be used:
**TODO:** add example.

Each aspect of `MaterialReconstruction.jl` is described in a dedicated section.

## Cost functions

Cost functions are the functions which are minimized during simulated
annealing. Usually, we have two binary arrays, one of them is fixed (the target
array) and the other is changed during annealing. An annealing process stops
when the value of the cost function becomes small enough which means two arrays
are "similar" to each other. A cost function takes two `CorrelationTracker`
objects (which are our arrays + a set of correlation functions calculated for
those arrays) and returns a real number which specifies how similar those
arrays are. Usually, this similarity is measured by comparing the values of
correlation functions. The following cost functions are defined in this package:

```@docs
euclid_mean
euclid_directional
euclid_mean_weighted
euclid_directional_weighted
čapek_cost
generalized_čapek_cost
```

## Initialization
