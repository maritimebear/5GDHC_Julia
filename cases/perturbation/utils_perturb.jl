# Utility functions to study time-evolution of perturbed states

module Utils_Perturb

export Polynomial_Tuple, PiecewisePolynomial, evalpiecewise


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


end # module
