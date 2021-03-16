
export BEMQuadAdaptive,
    BEMQuadQBF_Adaptive,
    BEMQuadQBF_Graded

using CompactTranslatesDict:
    PeriodicInterval

using BasisFunctions:
    unsafe_eval_element,
    unsafe_weight


using DomainIntegrals:
    QuadratureStrategy,
    QuadAdaptive,
    Q_quadgk,
    Q_hcubature,
    UnitIntervalRule,
    graded_rule_left,
    graded_rule_right,
    PointSingularity


###########################################
# Definition of BEM quadrature strategies
###########################################

"""
A `BEMQuadratureStrategy` groups quadrature strategies (as used by `DomainIntegrals.jl`)
for the regular, singular and nearly singular integrals corresponding to BEM matrix entries.

A concrete type should implement `regular`, `singular` and `nearlysingular` to
return the corresponding quadrature strategies.
"""
abstract type BEMQuadratureStrategy <: DomainIntegrals.QuadratureStrategy end

# For DomainIntegral strategies, there is no difference between the three types
regular(quad::QuadratureStrategy) = quad
singular(quad::QuadratureStrategy) = quad
nearlysingular(quad::QuadratureStrategy) = quad



"Generic adaptive quadrature for all BEM matrix entries."
struct BEMQuadAdaptive{T} <: BEMQuadratureStrategy
    adaptive    ::  QuadAdaptive{T}
end

regular(quad::BEMQuadAdaptive) = quad.adaptive
singular(quad::BEMQuadAdaptive) = quad.adaptive
nearlysingular(quad::BEMQuadAdaptive) = quad.adaptive

# Accept the same arguments as the QuadAdaptive constructor of DomainIntegrals.jl
BEMQuadAdaptive(args...) = BEMQuadAdaptive(QuadAdaptive(args...))
BEMQuadAdaptive{T}(atol::Number, args...) where {T} = BEMQuadAdaptive{T}(QuadAdaptive{T}(atol, args...))


"Combination of QBF and adaptive quadrature."
struct BEMQuadQBF_Adaptive{T} <: BEMQuadratureStrategy
    regular     ::  QuadQBF{T}
    adaptive    ::  QuadAdaptive{T}
end

BEMQuadQBF_Adaptive(regular::QuadQBF{T}, args...) where {T} =
    BEMQuadQBF_Adaptive(regular, QuadAdaptive{T}(args...))

regular(quad::BEMQuadQBF_Adaptive) = quad.regular
singular(quad::BEMQuadQBF_Adaptive) = quad.adaptive
nearlysingular(quad::BEMQuadQBF_Adaptive) = quad.adaptive

splineorder(q::BEMQuadQBF_Adaptive) = splineorder(q.regular)


"Combination of QBF and graded quadrature."
struct BEMQuadQBF_Graded{T} <: BEMQuadratureStrategy
    regular         ::  QuadQBF{T}
    graded_left     ::  UnitIntervalRule{T}
    graded_right    ::  UnitIntervalRule{T}
end

BEMQuadQBF_Graded(regular::QuadQBF{T}, σ = T(3)/10, n = 1, μ = one(T)/2, args...) where {T} =
    BEMQuadQBF_Graded(regular, graded_rule_left(σ, n, μ, args...), graded_rule_right(σ, n, μ, args...))

regular(quad::BEMQuadQBF_Graded) = quad.regular
singular(quad::BEMQuadQBF_Graded) = quad
nearlysingular(quad::BEMQuadQBF_Graded) = quad

splineorder(q::BEMQuadQBF_Graded) = splineorder(q.regular)


"A generic combination of BEM quadrature strategies."
struct BEMQuad <: BEMQuadratureStrategy
    regular
    singular
    nearlysingular
end

regular(quad::BEMQuad) = quad.regular
singular(quad::BEMQuad) = quad.singular
nearlysingular(quad::BEMQuad) = quad.nearlysingular


########################################################################
# Definition of projectionintegral and doubleprojectionintegral
########################################################################


function projectionintegral(qs, f, dict, idx, measure, sing = NoSingularity(), domain = support(dict, idx))
    integrand = t -> f(t)*unsafe_eval_element(dict, idx, t)*unsafe_weight(measure, t)
    integral(qs, integrand, domain, sing)
end

function projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing, domain)
    # For QBF we don't include the basis function
    integrand = t -> f(t)*unsafe_weight(measure, t)
    integral(qs, integrand, domain, sing)*sqrt(length(dict))
end

# Teach DomainIntegrals how to evaluate on a PeriodicInterval
function DomainIntegrals.integrate_domain(qs, integrand, domain::PeriodicInterval, measure, sing)
    if numelements(domain) > 1
        IEs = DomainIntegrals.integrate.(Ref(qs), Ref(integrand), elements(domain), Ref(measure), Ref(sing))
        DomainIntegrals.recombine_outcome(IEs)
    else
        DomainIntegrals.integrate_domain(qs, integrand, element(domain,1), measure, sing)
    end
end

# Special treatment for QBF quadrature: don't split the domain, instead periodize the integrand
function DomainIntegrals.integrate_domain(qs::QuadQBF, integrand, domain::PeriodicInterval, measure, sing)
    if numelements(domain) == 1
        DomainIntegrals.integrate_domain(qs, integrand, element(domain, 1), measure, sing)
    else
        # only periodize when crossing the boundary
        A, B = extrema(domain.periodicdomain)
        DomainIntegrals.integrate_domain(qs, t -> integrand(A+mod(t-A,B-A)), domain.subdomain, measure, sing)
    end
end

# Teach DomainIntegrals how to split a PeriodicInterval
function DomainIntegrals.splitdomain_point(x, domain::PeriodicInterval)
    if numelements(domain) == 1
        DomainIntegrals.splitdomain_point(x, element(domain, 1))
    else
        vcat(DomainIntegrals.splitdomain_point(x, element(domain, 1))..., DomainIntegrals.splitdomain_point(x, element(domain, 2))...)
    end
end


"Split a domain (originating from a spline) in subsets on which the spline is smooth."
piecewisedomains(domain::MappedDomain, order) = map_domain.(Ref(inverse_map(domain)), piecewisedomains(superdomain(domain), order))
function piecewisedomains(domain::AbstractInterval, order)
    a, b = extrema(domain)
    h = (b-a)/(order+1)
    [a+i*h..a+(i+1)*h for i in 0:order]
end
function piecewisedomains(domain::PeriodicInterval, order)
    realdomain = domain.subdomain
    periodicdomain = domain.periodicdomain
    a, b = extrema(realdomain)
    h = (b-a)/(order+1)
    segments = [a+i*h..a+(i+1)*h for i in 0:order]
    vcat([collect(elements(PeriodicInterval(s, periodicdomain))) for s in segments]...)
end

function projectionintegral(qs::BEMQuadQBF_Graded, f, dict, idx, measure, sing, domain)
    domains = piecewisedomains(support(dict, idx), splineorder(qs))
    sum(projectionintegral(qs, f, dict, idx, measure, sing, d) for d in domains)
end

function leftofsingularity(x, domain, support)
    L = width(support)
    if x - leftendpoint(domain) > L/2
        x = x - L
    elseif rightendpoint(domain) - x > L/2
        x = x + L
    end
    leftendpoint(domain) < x
end


function projectionintegral(qs::BEMQuadQBF_Graded, f, dict, idx, measure, sing::PointSingularity, domain::AbstractInterval)
    x = DomainIntegrals.point(sing)
    domains = DomainIntegrals.splitdomain_point(x, domain)
    z = zero(typeof(f(x)))
    for d in domains
        if leftofsingularity(x, d, support(dict))
            z += projectionintegral(qs.graded_right, f, dict, idx, measure, sing, d)
        else
            z += projectionintegral(qs.graded_left, f, dict, idx, measure, sing, d)
        end
    end
    z
end


doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing = NoSingularity()) =
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, support(dict1, idx1), support(dict2, idx2))

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
            domain1::AbstractInterval, domain2::AbstractInterval)
    integrand = t -> f(t[2],t[1])*unsafe_eval_element(dict1, idx1, t[1])*conj(unsafe_eval_element(dict2,idx2, t[2]))*
        unsafe_weight(measure1, t[1])*unsafe_weight(measure2, t[2])
    integral(qs, integrand, domain1 × domain2, sing)
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicInterval, domain2::PeriodicInterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, d1, d2)
            for d1 in elements(domain1), d2 in elements(domain2))
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicInterval, domain2::AbstractInterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, d1, domain2)
            for d1 in elements(domain1))
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::AbstractInterval, domain2::PeriodicInterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, domain1, d2)
            for d2 in elements(domain2))
end



doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
            domain1::AbstractInterval, domain2::AbstractInterval) =
        qbf_quadrature2(f, dict1, idx1, measure1, domain1, dict2, idx2, measure2, domain2,
            leftpoint(qs), rightpoint(qs), points(qs), weights(qs))

function doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicInterval, domain2::PeriodicInterval)

    A1, B1 = extrema(domain1.periodicdomain)
    A2, B2 = extrema(domain2.periodicdomain)
    f2 = t -> f(A1+mod(t[1]-A1,B1-A1), A2+mod(t[2]-A2,B2-A2))
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, domain1.subdomain, domain2.subdomain)
end
