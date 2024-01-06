module DHG

include("./FVM.jl")                     # provides submodule FVM
include("./TransportProperties.jl")     # submodule TransportProperties
include("./Utilities.jl")               # submodule Utilities
include("./DynamicalFunctions.jl")      # submodule DynamicalFunctions
include("./WrapperFunctions.jl")        # submodule WrapperFunctions
include("./NetworkComponents.jl")       # submodule NetworkComponents
include("./GraphParsing.jl")            # submodule GraphParsing
include("./ParameterStructs.jl")        # submodule ParameterStructs

# TODO: export parse_gml, others?
# TODO: global size_t, but NetworkDynamics takes dims::Int

end # module DHG
