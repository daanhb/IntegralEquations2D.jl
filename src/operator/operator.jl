
"""
The supertype of all integral operators.

An integral operator is characterized by the way it acts on a function:
integration with a kernel function and a certain integration measure.
"""
abstract type IntegralOperator <: AbstractOperator end

BasisFunctions.name(op::IntegralOperator) = "Integral operator"


# kernels typically have diagonal singularities.
struct LogDiagonalSingularity <: CurveSingularity end
struct RadialDiagonalSingularity <: CurveSingularity end
struct UnknownDiagonalSingularity <: CurveSingularity end

"Singularity of the kernel of the integral operator."
singularity(op::IntegralOperator) = UnknownDiagonalSingularity()



"The combination of an integral operator with a sampling operator."
struct SampledIntegralOperator <: AbstractOperator
    intop       ::  IntegralOperator
    sampling    ::  AbstractOperator
end

samplingoperator(op::SampledIntegralOperator) = op.sampling
integraloperator(op::SampledIntegralOperator) = op.intop

Base.:*(op::GridSampling, intop::IntegralOperator) = SampledIntegralOperator(intop, op)
Base.:*(op::ProjectionSampling, intop::IntegralOperator) = SampledIntegralOperator(intop, op)

BasisFunctions.name(op::SampledIntegralOperator) = _name(op, op.intop, op.sampling)

_name(::SampledIntegralOperator, intop, ::GridSampling) = "Collocation integral operator"
_name(::SampledIntegralOperator, intop, ::ProjectionSampling) = "Galerkin integral operator"


"The lazy application of an integral operator applied to a function."
struct LazyIntegralOperator
    intop   ::  IntegralOperator
    f
end

Base.:*(op::IntegralOperator, f) = LazyIntegralOperator(op, f)

(op::LazyIntegralOperator)(x) = eval_integral_operator(op.intop, op.f, x)

eval_integral_operator(op::IntegralOperator, f, x; options...) =
    _eval_integral_operator(kernel(op), domain(op), f, x; options...)

_eval_integral_operator(k::Kernel, domain, f, x; options...) =
    integral(t->eval_kernel(k, x, t) * f(t), domain; options...)
