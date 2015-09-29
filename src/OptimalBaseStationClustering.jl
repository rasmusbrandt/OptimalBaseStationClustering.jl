module OptimalBaseStationClustering

using CoordinatedPrecoding
import Lumberjack

export
    # assignment
    GeneralBranchAndBoundClustering,
    GeneralGreedyClustering,

    # precoding
    NoPrecoding

include("assignment/assignment.jl")
include("assignment/GeneralBranchAndBoundClustering.jl")
include("assignment/GeneralGreedyClustering.jl")
include("precoding/NoPrecoding.jl")

end # module
