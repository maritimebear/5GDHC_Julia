module Discretisation

export DiscretisationScheme, FVM
export upwind


abstract type DiscretisationScheme
    # Expected interface:
    # dx::Float64
    # convection::Function1
end


Base.@kwdef struct FVM{T1 <: Function} <: DiscretisationScheme
    dx::Float64         # discretisation sizing for temperature transport [m]
    convection::T1      # function implementing convection term in temperature transport equation
    # diffusion::T2
end


function upwind(phi::AbstractVector, phi_W, phi_E, u)
    # First-order upwind scheme
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


function linear(phi::AbstractVector, phi_W, phi_E, u)
    # Central difference scheme
    # Calculates closed surface integral (u*phi) . dS
    # Expects phi::Vector, where each element contains the value of phi in a finite-volume cell
    # Returns vector of results for each cell
    result = similar(phi)
    result[1] = 0.5 * (phi[1] + phi[2]) - phi_W # West boundary, Dirichlet BC: phi(west face) = phi_W
    result[end] = 0.5 * (phi[end] - phi[end-1]) # East boundary, Neumann BC: dphi/dx = 0 at east face
    for i in 2:(length(phi) - 1) # Interior
        result[i] = 0.5 * (phi[i+1] - phi[i-1])
    end
    return u .* result
end


function linear_upwind(phi::AbstractVector, phi_W, phi_E, u)
    # Linear upwind / second-order upwind scheme
    # Calculates closed surface integral (u*phi) . dS
    # Expects phi::Vector, where each element contains the value of phi in a finite-volume cell
    # Returns vector of results for each cell
    result = similar(phi)
    if u > 0
        result[1] = 1.5 * (phi[1] - phi_W)
        result[2] = (1.5 * phi[2]) - (2.0 * phi[1]) + (0.5 * phi_W)
        for i in 3:length(phi)
            result[i] = (1.5 * phi[i]) - (2.0 * phi[i-1]) + (0.5 * phi[i-2])
        end

    elseif u < 0
        result[end] = 1.5 * (phi_E - phi[end])
        result[end-1] = (2.0 * phi[end]) - (1.5 * phi[end-1]) - (0.5 * phi_E)
        for i in 1:(length(phi) - 2)
            result[i] = (2.0 * phi[i+1]) - (1.5 * phi[i]) - (0.5 * phi[i+2])
        end
    end
    return u .* result
end


end # module
