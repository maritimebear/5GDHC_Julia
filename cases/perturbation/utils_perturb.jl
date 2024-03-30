# Utility functions to study time-evolution of perturbed states

module Utils_Perturb

export Polynomial_Tuple, PiecewisePolynomial, evalpiecewise
export polynomial_coeffs, breakpoints, piecewise_linear


# Piecewise polynomials for prosumer controls

## "Performant piecewise function evaluation"
## https://discourse.julialang.org/t/performant-piecewise-function-evaluation/45762/3
## Answer by 'lucas711642'

struct Polynomial_Tuple{N, T <: Real}
	p::NTuple{N, T} # Polynomial coefficients: (c0, c1, c2, ...) for y = c0 + c1*x + c2*x^2 + ...
end


function (poly::Polynomial_Tuple)(x)
    return evalpoly(x, poly.p)
end


struct PiecewisePolynomial{N, P<:Polynomial_Tuple, T<:Real} <: Function
    polynomials::Vector{P} # Polynomial coefficients for each piece: [(c0, c1, c2, ...), ...]
    breakpoints::Vector{T} # Intersection points between each piece
	function PiecewisePolynomial{N}(polynomials::Vector{P}, breakpoints::Vector{T}) where {N, P <: Polynomial_Tuple, T <: Real}
		@assert length(polynomials) == N + 1 # Retain @assert here, let compiler decide whether to check
		@assert length(breakpoints) == N
		@assert issorted(breakpoints)
		new{N, P, T}(polynomials, breakpoints)
	end
end


function (p::PiecewisePolynomial{N})(x) where N
    return evalpiecewise(Val(N), p.polynomials, p.breakpoints, x)
end


function evalpiecewise(::Val{N}, functions, breakpoints, x) where N
	if @generated
		generator = (:(if x < breakpoints[$k]
			return functions[$k](x)
		end) for k = 1 : N)
		quote
			@inbounds begin
				$(generator...)
				return functions[$(N + 1)](x)
			end
		end
	else
		ind = searchsortedfirst(breakpoints, x)
		return functions[ind](x)
	end
end


# Convenience functions to construct piecewise polynomials

CoordinatesTuple = NTuple{2, <:Real} # type alias for convenience, tuple of (x, y) coordinates


function breakpoints(xys::CoordinatesTuple...)
    # -> Vector{Float64}
    # Intended to be called by other functions to define piecewise polynomials
    # Returns x-coordinates of points where pieces intersect

    return [Float64(xy[1]) for xy in xys[2:end-1]]
        # First and last (x, y) are not intersection points between pieces
end


function breakpoints(_::Vararg{CoordinatesTuple, 2})
    # -> Nothing
    # Method to catch edge case, throws exception
    error("Number of (x, y) pairs must be > 2")
end


function breakpoints(_::CoordinatesTuple)
    # -> Nothing
    # Method to catch edge case, throws exception
    error("Number of (x, y) pairs must be > 2")
end


# Wrapper functions to construct piecewise-linear polynomials for prosumer controls

function polynomial_coeffs(xy1::CoordinatesTuple, xy2::CoordinatesTuple)
    # -> Polynomial_Tuple{2, Float64}: (c0, c1) for y = c0 + (c1 * x)
    # Calculates polynomial coefficients for function passing through each specified xy-coordinate pair
    # Input type: (x, y), each xy pair specifies coordinates of fit points
    # Method specialisation for linear polynomials, implemented by type signature (2 xy pairs => linear polynomial)

    x1 = xy1[1]; x2 = xy2[1]; y1 = xy1[2]; y2 = xy2[2];
    c1 = (y2 - y1) / (x2 - x1)
    c0 = y1 - (x1 * c1)
    return Polynomial_Tuple((c0, c1))
end


function piecewise_linear(xys::Vararg{CoordinatesTuple, N}) where {N}
    # -> PiecewisePolynomial{N-2, PolynomialTuple{2, Float64}, Float64}
    # Constructs piecewise-linear polynomial between specified (x, y) points
    # Input type: (x, y), each xy pair specifies fit points for polynomial
    # First and last (x, y) must be boundary points

    coefficients = [polynomial_coeffs(xys[i], xys[i+1]) for i in 1:(N - 1)]
    intersections = breakpoints(xys...) # Will throw if N < 3
    return PiecewisePolynomial{N-2}(coefficients, intersections)
end


end # module
