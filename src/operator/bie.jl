
"A measure associated with a parameterization "
struct ParameterizationMeasure{P,T} <: DomainIntegrals.Weight{T}
    param   ::  P
end

ParameterizationMeasure(param) = ParameterizationMeasure{typeof(param),eltype(domain(param))}(param)

DomainIntegrals.unsafe_weightfun(p::ParameterizationMeasure, t) = gradientnorm(p.param, t)
DomainIntegrals.support(p::ParameterizationMeasure) = domain(p.param)


export BoundaryIntegralOperator
"A boundary integral operator with a parameterized boundary."
struct BoundaryIntegralOperator{K,P,T} <: IntegralOperator
    kernel      ::  K
    boundary    ::  Domain{T}
    param       ::  P
end

BasisFunctions.name(op::BoundaryIntegralOperator) = "Boundary integral operator"

kernel(op::BoundaryIntegralOperator) = op.kernel
boundary(op::BoundaryIntegralOperator) = op.boundary
parameterization(op::BoundaryIntegralOperator) = op.param
parameterdomain(op::BoundaryIntegralOperator) = domain(parameterization(op))

eval_kernel(op::BoundaryIntegralOperator, x, y) =
    kernel(op)(x, y, parameterization(op))

measure(op::BoundaryIntegralOperator) = ParameterizationMeasure(parameterization(op))

singularity(op::BoundaryIntegralOperator) = singularity(kernel(op))

export op_apply
"Apply the integral operator to the given density function and evaluate at the point t."
function op_apply(op::BoundaryIntegralOperator, density, t)
    μ = measure(op)
    integrand = y -> eval_kernel(op, t, y) * density(y) * unsafe_weightfun(μ, y)
    DomainIntegrals.integral(QuadAdaptive(), integrand, parameterdomain(op), LogSingPoint(t))
end
