using Documenter
using CorrelationFunctions
using CorrelationTrackers
using MaterialReconstruction

makedocs(sitename = "MaterialReconstruction.jl documentation",
         format = Documenter.HTML(prettyurls = false))

deploydocs(repo = "github.com/fatimp/MaterialReconstruction.jl.git")
