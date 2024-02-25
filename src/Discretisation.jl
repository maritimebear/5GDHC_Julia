module Discretisation

export DiscretisationScheme, FVM
export upwind


abstract type DiscretisationScheme end

Base.@kwdef struct FVM{T1 <: Function} <: DiscretisationScheme
    dx::Float64         # discretisation sizing for temperature transport [m]
    convection::T1      # function implementing convection term in temperature transport equation
    # diffusion::T2
end


function upwind(phi::AbstractVector, phi_W, phi_E, u)
    # Calculates closed surface integral (u*phi) . dS
    # Expects phi::Vector, where each element contains the value of phi in a finite-volume cell
    # Returns vector of results for each cell

    neighbour = similar(phi)
    if u > 0
        neighbour[1] = phi_W
        neighbour[2:end] .= phi[1:end-1]
    elseif u < 0
        neighbour[1:end-1] .= phi[2:end]
        neighbour[end] = phi_E
    end
    return (abs(u) .* (phi .- neighbour))
end

end # module
