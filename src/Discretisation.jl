module Discretisation

export DiscretisationScheme, FVM

export upwind, linear, linear_upwind
export create_TVD_scheme
export vanLeer, vanAlbada, minmod, lu_TVD


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


# TVD schemes

function grad_ratio(phi::AbstractVector, phi_W, phi_E, u)
    # -> Vector{Float64}
    # Calculates gradient ratio required by flux limiters in TVD schemes
    result = similar(phi, length(phi)+1) # length(result) == number of faces in grid
    if u >= 0
        result[1] = 0.0 # Interface with west-boundary node
        result[2] = (phi[1] - phi_W) / (phi[2] - phi[1])
        for i in 3:(length(result) - 1)
            result[i] = (phi[i-1] - phi[i-2]) / (phi[i] - phi[i-1])
        end
        result[end] = (phi[end] - phi[end-1]) / (phi_E - phi[end]) # Interface with east-boundary node

    else # u < 0
        result[1] = (phi[1] - phi[2]) / (phi_W - phi[1]) # Interface with west-boundary node
        for i in 2:(length(result) - 2)
            result[i] = (phi[i] - phi[i+1]) / (phi[i-1] - phi[i])
        end
        result[end-1] = (phi[end] - phi_E) / (phi[end-1] - phi[end])
        result[end] = 0.0 # Interface with east-boundary node
    end
    # At locations where phi is constant, the gradient ratio = NaN or Inf or -Inf due to division by zero
    # These NaNs are replaced with zeros, as at these locations the flux limiter can also go to zero (i.e. Upwind scheme)
    # This does not cause a loss in accuracy, as the function phi is constant anyway [https://folk.ntnu.no/leifh/teaching/tkt4140/._main074.html]
    map!(x -> (isfinite(x) ? x : 0.0), result, result)
    return result
end


function create_TVD_scheme(limiter_function)
    # Returns closure that implements a TVD convection scheme with the specified flux limiter function
    function convection_scheme(phi::AbstractVector, phi_W, phi_E, u)
        # Closure, flux limiter captured from outer function
        let limiter_fn = limiter_function
            grad_ratios = grad_ratio(phi, phi_W, phi_E, u)
            phi_upstream = u >= 0 ? [phi_W; phi[1:end]] : [phi[1:end]; phi_E]
            phi_downstream = u >= 0 ? [phi[1:end]; phi_E] : [phi_W; phi[1:end]]
            face_phis = phi_upstream .+ (0.5 .* limiter_fn.(grad_ratios) .* (phi_downstream .- phi_upstream))
            return u .* (face_phis[2:end] .- face_phis[1:end-1])
        end # let block
    end # Closure
    return convection_scheme
end


function fluxlimiter_vanLeer(grad_ratio)
    # -> Vector{Float64}
    # van Leer flux limiter, calculates limit for each face in grid
    return (grad_ratio .+ abs.(grad_ratio)) ./ (1.0 .+ abs.(grad_ratio))
end


## TVD schemes as closures

vanLeer = create_TVD_scheme((r) -> ((r .+ abs.(r)) ./ (1.0 .+ abs.(r))))

vanAlbada = create_TVD_scheme((r) -> ((r.^2 .+ r) ./ (r.^2 .+ 1.0)))

minmod = create_TVD_scheme((r) -> (clamp.(r, 0.0, 1.0)))

lu_TVD = create_TVD_scheme((_) -> 1.0) # TVD version of linear upwind, for testing purposes


end # module
