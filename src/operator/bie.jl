
"An integral operator defined by a kernel function of two variables."
struct KernelIntegralOperator{S,T,K<:Kernel{S,T}} <: IntegralOperator{S,T}
    kernel  ::  K
    domain  ::  Domain{S}
end

kernel(op::KernelIntegralOperator) = op.kernel
domain(op::KernelIntegralOperator) = op.domain

eval_kernel(op::KernelIntegralOperator, x, y) = eval_kernel(op.kernel, x, y)

measure(op::KernelIntegralOperator{S}) where {S} = BasisFunctions.WholeLebesgueMeasure{S}()

singularity(op::KernelIntegralOperator) = singularity(kernel(op))


struct ParameterizationMeasure{P,S,T} <: BasisFunctions.Measure{S}
    param   ::  P
end
BasisFunctions.rangetype(p::ParameterizationMeasure{P,S,T}, x) where {P,S,T} = T
BasisFunctions.unsafe_weight(p::ParameterizationMeasure, x) = gradientnorm(p.param, x)
BasisFunctions.support(p::ParameterizationMeasure) = domain(p.param)


"A boundary integral operator with a parameterized boundary."
struct BoundaryIntegralOperator{S,T,P,K<:Kernel{S,T}} <: IntegralOperator{S,T}
    kernel      ::  K
    domain      ::  Domain{S}
    param       ::  P
    paramdomain ::  Domain
end

BoundaryIntegralOperator(kernel, k_domain, param) =
    BoundaryIntegralOperator(kernel, k_domain, param, domain(param))

BoundaryIntegralOperator(kernel::Kernel{S,T}, domain, param, paramdomain) where {S,T} =
    BoundaryIntegralOperator{S,T,typeof(param),typeof(kernel)}(kernel, domain, param, paramdomain)

BasisFunctions.name(op::BoundaryIntegralOperator) = "Boundary integral operator"

kernel(op::BoundaryIntegralOperator) = op.kernel
domain(op::BoundaryIntegralOperator) = op.domain
parameterization(op::BoundaryIntegralOperator) = op.param
paramdomain(op::BoundaryIntegralOperator) = op.paramdomain

eval_kernel(op::BoundaryIntegralOperator, t, tau) =
    _eval_kernel(op, t, tau, op.param, op.kernel)

eval_kernel(op::BoundaryIntegralOperator, t, tau, x, y) =
    _eval_kernel(op, t, tau, op.param, op.kernel, x, y)

_eval_kernel(op::BoundaryIntegralOperator, t, tau, param, kernel,
    x = applymap(param, t), y = applymap(param, tau)) = eval_kernel(kernel, x, y)

measure(op::BoundaryIntegralOperator{S,T,P}) where {S,T,P} = ParameterizationMeasure{P,S,T}(op.param)

isperiodic(op::BoundaryIntegralOperator) = paramdomain(op) isa Interval{:closed,:open}

period(op::BoundaryIntegralOperator) = width(paramdomain(op))

singularity(op::BoundaryIntegralOperator) = singularity(kernel(op))


Base.:*(op::GridSampling, iop::IntegralOperator) = CollocationIntOp(op, iop)
Base.:*(op::ProjectionSampling, iop::IntegralOperator) = GalerkinIntOp(op, iop)
