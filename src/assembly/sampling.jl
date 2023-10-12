
"""
A `QuadProjectionSampling` is an operator that maps a function to its inner products
with a projection basis, with those inner products computed by a quadrature strategy.
"""
struct QuadProjectionSampling <: BasisFunctions.SamplingOperator
    dict		::  Dictionary
	src_space	::	BasisFunctions.FunctionSpace
	quad		::	QuadratureStrategy
end

QuadProjectionSampling(dict::Dictionary, quad::QuadratureStrategy) = QuadProjectionSampling(dict, Span(dict), quad)

# BasisFunctions.name(op::QuadProjectionSampling) = "Projection by quadrature operator"

BasisFunctions.dictionary(op::QuadProjectionSampling) = op.dict

BasisFunctions.dest(op::QuadProjectionSampling) = dictionary(op)

BasisFunctions.src_space(op::QuadProjectionSampling) = op.src_space

quad(op::QuadProjectionSampling) = op.quad

BasisFunctions.apply!(result, op::QuadProjectionSampling, f) =
	quadproject!(result, f, dictionary(op), quad(op))

function quadproject!(result, f, dict::Dictionary1d, quad::QuadratureStrategy)
	@assert size(result) == size(dict)
    for i in eachindex(dict)
		result[i] = projectionintegral(quad, f, dict, i)
	end
	result
end
