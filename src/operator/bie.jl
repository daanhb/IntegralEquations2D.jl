
"A measure associated with a parameterization "
struct JacobianMeasure{P,T} <: BasisFunctions.Measure{T}
    param   ::  P
end

JacobianMeasure(param) = JacobianMeasure{typeof(param),eltype(domain(param))}(param)

BasisFunctions.unsafe_weight(p::JacobianMeasure, x) = jacobian(p.param, x)
BasisFunctions.support(p::JacobianMeasure) = domain(p.param)


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

measure(op::BoundaryIntegralOperator) = JacobianMeasure(parameterization(op))

singularity(op::BoundaryIntegralOperator) = singularity(kernel(op))
