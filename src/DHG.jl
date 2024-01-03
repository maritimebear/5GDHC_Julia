module DHG

import NetworkDynamics
import LinearAlgebra

include("./InterpolationSchemes.jl")    # provides submodule FVM
include("./Utilities.jl")               # submodule Utilities
include("./Physics.jl")                 # submodule Physics
include("./WrapperFunctions.jl")        # submodule WrapperFunctions
include("./NetworkComponents.jl")       # submodule NetworkComponents
include("./GraphParsing.jl")            # submodule GraphParsing
include("./ParameterStructs.jl")        # submodule ParameterStructs

end # module DHG
