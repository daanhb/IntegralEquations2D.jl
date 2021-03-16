
export eval_field
"Evaluate the field at a point x"
eval_field(A::BEMOperator, density::Expansion, x, quad = A.quad) =
    eval_field(integraloperator(A), density, x, quad)

function eval_field(intop::BoundaryIntegralOperator, density::Expansion, x, quad)
    integrand = t -> eval_kernel(intop, x, t) * unsafe_weightfun(measure(intop), t) * density(t)
    DomainIntegrals.integral(quad, integrand, support(density))
end

eval_field(intop::BoundaryIntegralOperator, density::Expansion, x, quad::QuadQBF) =
    _eval_field(intop, dictionary(density), coefficients(density), x, quad)

function _eval_field(intop, dict, coef, x, quad::QuadQBF)
    z = zero(eltype(coef))
    for i in 1:length(coef)
        z = z + coef[i] * collocation_BEM_entry(x, (dict, i), intop, quad)
    end
    z
end
