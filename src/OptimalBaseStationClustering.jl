VERSION >= v"0.4.0-dev+6521" && __precompile__()

module OptimalBaseStationClustering

using CoordinatedPrecoding
import Lumberjack, PyPlot

export
    # assignment
    GeneralBranchAndBoundClustering,
    GeneralGreedyClustering,

    # misc
    display_clustering

include("assignment/assignment.jl")
include("assignment/GeneralBranchAndBoundClustering.jl")
include("assignment/GeneralGreedyClustering.jl")
include("misc/visualization.jl")

end # module
