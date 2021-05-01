
export compute_BEM_entry

# In this file we implement the `compute_BEM_entry` function.
# This function determines whether the corresponding integral is regular,
# singular or nearly singular. It differentiates between collocation and
# Galerkin discretization to do so.
# In the end, the function invokes `projectionintegral` (for collocation)
# or `doubleprojectionintegral` (for Galerkin).

"Compute the BEM matrix entry A[i,j]."
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


"Does the collocation entry correspond to a (nearly) singular integral, up to the given distance threshold `δ`?"
issingular_coll(sing, x, ϕ_j) =
    issingular_coll(sing, x, ϕ_j, DomainSets.default_tolerance(support(ϕ_j...)))

issingular_coll(::NoSingularity, x, ϕ_j, δ) = false

issingular_coll(sing, x, ϕ_j, δ) =
    _issingular_coll(x, support(ϕ_j...), δ)
_issingular_coll(x, d::AbstractInterval, δ) = approx_in(x, d, δ)
_issingular_coll(x, d::MappedDomain, δ) = _issingular_coll(inverse_map(d, x), superdomain(d), δ)
# Line below catches a corner case
_issingular_coll(x, d::PeriodicInterval, δ) =
    approx_in(x, d, δ) || approx_in(x+period(d), d, δ) || approx_in(x-period(d), d, δ)

"Which threshold decides whether an integral is nearly singular?"
nearsingular_threshold(ϕ_j) = nearsingular_threshold(support(ϕ_j...))
# By default we choose twice the length of a single basis function
nearsingular_threshold(d::AbstractInterval) = 2DomainSets.width(d)
nearsingular_threshold(d::PeriodicInterval) = 2sum(map(DomainSets.width, components(d)))
# we don't know the radius of a general domain...
nearsingular_threshold(d::Domain) = 0.01

function collocation_BEM_entry(x, ϕ_j, intop, quad, sing = singularity(intop))
    if issingular_coll(sing, x, ϕ_j)
        singular_collocation_entry(singular(quad), x, ϕ_j, intop, sing)
    elseif issingular_coll(sing, x, ϕ_j, nearsingular_threshold(ϕ_j))
        nearly_singular_collocation_entry(nearlysingular(quad), x, ϕ_j, intop, sing)
    else
        regular_collocation_entry(regular(quad), x, ϕ_j, intop)
    end
end

collocation_singularity(x, sing) = sing
collocation_singularity(x, sing::LogSingularDiagonal) = LogSingPoint(x)

singular_collocation_entry(quad, x, ϕ_j, intop, sing) =
    projectionintegral(singular(quad), y -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), collocation_singularity(x, sing))

nearly_singular_collocation_entry(quad, x, ϕ_j, intop, sing) =
    projectionintegral(nearlysingular(quad), y -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), collocation_singularity(x, sing))

regular_collocation_entry(quad, x, ϕ_j, intop) =
    projectionintegral(regular(quad), y -> eval_kernel(intop, x, y), ϕ_j..., measure(intop))


##############
# Galerkin
##############

"Does the Galerkin entry correspond to a (nearly) singular integral, up to the given distance threshold `δ`?"
issingular_galerkin(sing, ϕ_i, ϕ_j) =
    issingular_galerkin(sing, ϕ_i, ϕ_j, DomainSets.default_tolerance(support(ϕ_j...)))

issingular_galerkin(::NoSingularity, ϕ_i, ϕ_j, δ) = false
issingular_galerkin(sing, ϕ_i, ϕ_j, δ) =
    _issingular_galerkin(sing, ϕ_i, ϕ_j, δ, support(ϕ_i...), support(ϕ_j...))
_issingular_galerkin(sing, ϕ_i, ϕ_j, δ, domain1, domain2) =
    distance(domain1, domain2) < δ

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
    min(distance(sub1, sub2), distance(sub1 .+ period,sub2), distance(sub1 .- period,sub2))
end

distance(domain1::UnitInterval, domain2::PeriodicInterval) = zero(eltype(domain1))

function galerkin_BEM_entry(ϕ_i, ϕ_j, μ_projection, intop, quad, sing = singularity(intop))
    if issingular_galerkin(sing, ϕ_i, ϕ_j)
        singular_galerkin_BEM_entry(singular(quad), sing, ϕ_i, μ_projection, ϕ_j, intop)
    elseif issingular_galerkin(sing, ϕ_i, ϕ_j, nearsingular_threshold(ϕ_j))
        nearly_singular_galerkin_BEM_entry(nearlysingular(quad), sing, ϕ_i, μ_projection, ϕ_j, intop)
    else
        regular_galerkin_BEM_entry(regular(quad), ϕ_i, μ_projection, ϕ_j, intop)
    end
end

galerkin_singularity(sing) = sing
galerkin_singularity(sing::LogSingularDiagonal) = SingularDiagonal()

singular_galerkin_BEM_entry(quad, sing, ϕ_i, μ_projection, ϕ_j, intop) =
    doubleprojection(singular(quad), (x,y) -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), ϕ_i..., μ_projection, galerkin_singularity(sing))

nearly_singular_galerkin_BEM_entry(quad, sing, ϕ_i, μ_projection, ϕ_j, intop) =
    doubleprojection(nearlysingular(quad), (x,y) -> eval_kernel(intop, x, y),
            ϕ_j..., measure(intop), ϕ_i..., μ_projection, galerkin_singularity(sing))

regular_galerkin_BEM_entry(quad, ϕ_i, μ_projection, ϕ_j, intop) =
    doubleprojection(regular(quad), (x,y) -> eval_kernel(intop, x, y), ϕ_j...,
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
