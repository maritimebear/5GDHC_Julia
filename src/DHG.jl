module DHG

include("Discretisation.jl")
include("Fluids.jl")
include("Transport.jl")
include("Miscellaneous.jl")
include("NetworkComponents.jl")
include("DynamicalFunctions.jl")
include("WrapperFunctions.jl")

import .NetworkComponents: Node, Edge, Prosumer, JunctionNode, ReferenceNode, Pipe, PressureChange, Massflow
import .Transport: TransportModels
import .WrapperFunctions: node, edge



end
