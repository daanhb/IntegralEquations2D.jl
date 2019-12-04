
using BasisFunctions:
    unsafe_eval_element

compute_BEM_entry(op::DenseBEMOperator, i::Int, j::Int, quad = op.quad) =
    compute_BEM_entry(i, (src(op), j), samplingoperator(op), integraloperator(op), quad)

# Catch collocation/galerkin based on the sampling operator
compute_BEM_entry(i, ϕ_j, sampling::GridSampling, intop, quad) =
    collocation_BEM_entry(grid(sampling), i, ϕ_j, intop, quad)

compute_BEM_entry(i, ϕ_j, sampling::ProjectionSampling, intop, quad) =
    galerkin_BEM_entry((dictionary(sampling), i), ϕ_j, measure(sampling), intop, quad)


##############
# Collocation
##############

collocation_BEM_entry(grid::AbstractGrid, i, ϕ_j, intop, quad) =
    collocation_BEM_entry(grid[i], ϕ_j, intop, quad)


# Does the collocation entry correspond to a singular integral?
issingular_coll(quad, ::NoSingularity, x, ϕ_j) = false
issingular_coll(quad, sing, x, ϕ_j) =
    _issingular_coll(x, support(ϕ_j...))
_issingular_coll(x, d::AbstractInterval) = approx_in(x, d)
# Line below catches a corner case
_issingular_coll(x, d::PeriodicInterval) =
    approx_in(x, d) || approx_in(x+period(d), d) || approx_in(x-period(d), d)

function collocation_BEM_entry(x, ϕ_j, intop, quad, sing = singularity(intop))
    if issingular_coll(quad, sing, x, ϕ_j)
        singular_collocation_entry(quad, x, ϕ_j, intop, sing)
    else
        regular_collocation_entry(quad, x, ϕ_j, intop)
    end
end

collocation_singularity(x, sing) = sing
collocation_singularity(x, sing::LogDiagonalSingularity) = LogPointSingularity(x)

singular_collocation_entry(qs, x, ϕ_j, intop, sing) =
    projectionintegral(qs, y -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), collocation_singularity(x, sing))

regular_collocation_entry(qs, x, ϕ_j, intop) =
    projectionintegral(qs, y -> eval_kernel(intop, x, y), ϕ_j..., measure(intop))


##############
# Galerkin
##############

issingular_galerkin(quad, ::NoSingularity, ϕ_i, ϕ_j) = false
issingular_galerkin(quad, sing, ϕ_i, ϕ_j) =
    _issingular_galerkin(quad, sing, ϕ_i, ϕ_j, support(ϕ_i...), support(ϕ_j...))
_issingular_galerkin(quad, sing, ϕ_i, ϕ_j, domain1, domain2) =
    distance(domain1, domain2) < 100eps(eltype(domain1))

function distance(domain1::AbstractInterval, domain2::AbstractInterval)
    a1 = leftendpoint(domain1)
    b1 = rightendpoint(domain1)
    a2 = leftendpoint(domain2)
    b2 = rightendpoint(domain2)
    if a2 > b1
        a2-b1
    elseif a1 > b2
        a1-b2
    else
        zero(a1)
    end
end

function distance(domain1::PeriodicInterval, domain2::PeriodicInterval)
    @assert domain1.periodicdomain == domain2.periodicdomain
    period = width(domain1.periodicdomain)
    sub1 = domain1.subdomain
    sub2 = domain2.subdomain
    min(distance(sub1, sub2), distance(sub1+period,sub2), distance(sub1-period,sub2))
end

function galerkin_BEM_entry(ϕ_i, ϕ_j, μ_projection, intop, quad, sing = singularity(intop))
    if issingular_galerkin(quad, sing, ϕ_i, ϕ_j)
        singular_galerkin_BEM_entry(quad, sing, ϕ_i, μ_projection, ϕ_j, intop)
    else
        regular_galerkin_BEM_entry(quad, ϕ_i, μ_projection, ϕ_j, intop)
    end
end

galerkin_singularity(sing) = sing
galerkin_singularity(sing::LogDiagonalSingularity) = DiagonalSingularity()

singular_galerkin_BEM_entry(quad, sing, ϕ_i, μ_projection, ϕ_j, intop) =
    doubleprojection(quad, (x,y) -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), ϕ_i..., μ_projection, galerkin_singularity(sing))

regular_galerkin_BEM_entry(quad, ϕ_i, μ_projection, ϕ_j, intop) =
    doubleprojection(quad, (x,y) -> eval_kernel(intop, x, y), ϕ_j...,
            measure(intop), ϕ_i..., μ_projection)


#######################
# Ranges and subblocks
#######################

# There is overhead in the computation of a single entry that can be
# saved if we are computing a range of entries.

compute_BEM_range!(result, A::DenseBEMOperator, I, J, quad = A.quad) =
    compute_BEM_range!(result, I, J, quad, A.src, A.sampling, A.intop)

compute_BEM_range!(result, I, J, quad, dict, sampling::GridSampling, intop) =
    collocation_BEM_range!(result, grid(sampling), I, quad, dict, J, intop, measure(intop))

collocation_BEM_range!(result, grid::AbstractGrid, i::Int, quad, dict, J, intop, measure) =
    collocation_BEM_range!(result, grid[i], i, quad, dict, J, intop, measure)

function collocation_BEM_range!(result, x, i, quad, dict, J, intop, measure)
    sing = singularity(intop)
    for j in J
        result[i,j] = collocation_BEM_entry(x, quad, dict, j, intop, measure, sing)
    end
    result
end

function collocation_BEM_range!(result, grid::AbstractGrid, I, quad, dict, J, intop, measure)
    for i in I
        collocation_BEM_range!(result, grid, i, quad, dict, J, intop, measure)
    end
    result
end
