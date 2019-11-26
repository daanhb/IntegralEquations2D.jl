
"An integral operator defined by a kernel function of two variables."
struct KernelIntOp{S,T,K<:Kernel{S,T}} <: IntegralOperator{S,T}
    kernel  ::  K
    domain  ::  Domain{S}
end

kernel(op::KernelIntOp) = op.kernel
domain(op::KernelIntOp) = op.domain

eval_kernel(op::KernelIntOp, x, y) = eval_kernel(op.kernel, x, y)

measure(op::KernelIntOp{S}) where {S} = BasisFunctions.WholeLebesgueMeasure{S}()

singularity(op::KernelIntOp) = singularity(kernel(op))


struct ParameterizationMeasure{P,S,T} <: BasisFunctions.Measure{S}
    param   ::  P
end
BasisFunctions.rangetype(p::ParameterizationMeasure{P,S,T}, x) where {P,S,T} = T
BasisFunctions.unsafe_weight(p::ParameterizationMeasure, x) = gradientnorm(p.param, x)
BasisFunctions.support(p::ParameterizationMeasure) = domain(p.param)


"A boundary integral operator with a parameterized boundary."
struct BoundaryIntOp{S,T,P,K<:Kernel{S,T}} <: IntegralOperator{S,T}
    kernel      ::  K
    domain      ::  Domain{S}
    param       ::  P
    paramdomain ::  Domain
end

BoundaryIntOp(kernel, k_domain, param) =
    BoundaryIntOp(kernel, k_domain, param, domain(param))

BoundaryIntOp(kernel::Kernel{S,T}, domain, param, paramdomain) where {S,T} =
    BoundaryIntOp{S,T,typeof(param),typeof(kernel)}(kernel, domain, param, paramdomain)

BasisFunctions.name(op::BoundaryIntOp) = "Boundary integral operator"

kernel(op::BoundaryIntOp) = op.kernel
domain(op::BoundaryIntOp) = op.domain
parameterization(op::BoundaryIntOp) = op.param
paramdomain(op::BoundaryIntOp) = op.paramdomain

eval_kernel(op::BoundaryIntOp, t, tau) =
    _eval_kernel(op, t, tau, op.param, op.kernel)

eval_kernel(op::BoundaryIntOp, t, tau, x, y) =
    _eval_kernel(op, t, tau, op.param, op.kernel, x, y)

_eval_kernel(op::BoundaryIntOp, t, tau, param, kernel,
    x = applymap(param, t), y = applymap(param, tau)) = eval_kernel(kernel, x, y)

measure(op::BoundaryIntOp{S,T,P}) where {S,T,P} = ParameterizationMeasure{P,S,T}(op.param)

isperiodic(op::BoundaryIntOp) = paramdomain(op) isa Interval{:closed,:open}

period(op::BoundaryIntOp) = width(paramdomain(op))

singularity(op::BoundaryIntOp) = singularity(kernel(op))


Base.:*(op::GridSampling, iop::IntegralOperator) = CollocationIntOp(op, iop)
Base.:*(op::ProjectionSampling, iop::IntegralOperator) = GalerkinIntOp(op, iop)
