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

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

end # module
