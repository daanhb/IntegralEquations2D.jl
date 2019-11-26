
"""
The supertype of all integral operators.

The kernel function maps arguments of type `S` to an output of type `T`.
"""
abstract type IntegralOperator{S,T} <: AbstractOperator end

BasisFunctions.name(op::IntegralOperator) = "Integral operator"



abstract type DiscretizedIntOp{S,T} <: AbstractOperator end

"Collocation integral operator"
struct CollocationIntOp{S,T} <: DiscretizedIntOp{S,T}
    sampling    ::  GridSampling
    intop       ::  IntegralOperator{S,T}
end

CollocationIntOp{S,T}(intop::IntegralOperator{S,T}, grid::AbstractGrid{S}) where {S,T} =
    CollocationIntOp{S,T}(intop, GridSampling(grid, T))

BasisFunctions.name(op::CollocationIntOp) = "Collocation integral operator"

"Galerkin integral operator"
struct GalerkinIntOp{S,T} <: DiscretizedIntOp{S,T}
    sampling    ::  ProjectionSampling
    intop       ::  IntegralOperator{S,T}
end

BasisFunctions.name(op::GalerkinIntOp) = "Galerkin integral operator"




"The lazy application of an integral operator applied to a function."
struct LazyIntOp
    intop   ::  IntegralOperator
    f
end

Base.:*(op::IntegralOperator, f) = LazyIntOp(op, f)

(op::LazyIntOp)(x) = eval_integral_operator(op.intop, op.f, x)

eval_integral_operator(op::IntegralOperator, f, x; options...) =
    _eval_integral_operator(kernel(op), domain(op), f, x; options...)

_eval_integral_operator(k::Kernel, domain, f, x; options...) =
    integral(t->eval_kernel(k, x, t) * f(t), domain; options...)
