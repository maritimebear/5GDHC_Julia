module DHG

include("./Utilities.jl")               # submodule Utilities
include("./FVM.jl")                     # provides submodule FVM
include("./NetworkComponents.jl")       # submodule NetworkComponents
include("./TransportProperties.jl")     # submodule TransportProperties
include("./DynamicalFunctions.jl")      # submodule DynamicalFunctions
include("./WrapperFunctions.jl")        # submodule WrapperFunctions
include("./GraphParsing.jl")            # submodule GraphParsing
include("./ParameterStructs.jl")        # submodule ParameterStructs

# TODO: export parse_gml, others?
# TODO: global size_t, but NetworkDynamics takes dims::Int

import .TransportProperties: TransportCoefficients
import .ParameterStructs: NodeParameters, EdgeParameters, GlobalParameters, Parameters

export TransportCoefficients
export NodeParameters, EdgeParameters, GlobalParameters, Parameters

end # module DHG
