
suggest_quadrature(A::BoundaryIntOp) = QuadAdaptive()

eval_field(A::BEMOperator, density::Expansion, x, quad = A.quad) =
    eval_field(integraloperator(A), density, x, quad)


function eval_field(intop::BoundaryIntOp, density::Expansion, x, quad = suggest_quadrature(intop))
    integrand = t -> eval_kernel(intop.kernel, x, applymap(intop.param, t)) * unsafe_weight(measure(intop), t) * density(t)
    DomainIntegrals.integral(integrand, support(density))
end
