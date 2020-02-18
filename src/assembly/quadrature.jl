
using CompactTranslatesDict:
    PeriodicInterval

using BasisFunctions:
    unsafe_eval_element,
    unsafe_weight

# # A few routines for interoperability with pretty printing
# BasisFunctions.hasstencil(q::QuadratureStrategy) = false
# BasisFunctions.symbol(q::QuadratureStrategy) = "Quad"
# BasisFunctions.iscomposite(q::QuadratureStrategy) = false

import DomainIntegrals:
    QuadratureStrategy, QuadAdaptive,
    Q_quadgk, Q_hcubature

projectionintegral(qs, f, dict, idx, measure, sing = NoSingularity()) =
    projectionintegral(qs, f, dict, idx, measure, sing, support(dict, idx))


function projectionintegral(qs, f, dict, idx, measure, sing, domain)
    if domain isa PeriodicInterval
        sum(projectionintegral(qs, f, dict, idx, measure, sing, el) for el in elements(domain))
    else
        integrand = t -> f(t)*unsafe_eval_element(dict, idx, t)*unsafe_weight(measure, t)
        DomainIntegrals.integral(qs, integrand, domain, sing)
    end
end

# # Teach DomainIntegrals how to evaluate on a PeriodicInterval
# function DomainIntegrals.quadrature_d(qs, integrand, domain::PeriodicInterval, measure, sing)
#     if numelements(domain) > 1
#         IEs = DomainIntegrals.quadrature.(Ref(qs), Ref(integrand), elements(domain), Ref(measure), Ref(sing))
#         DomainIntegrals.recombine_outcome(IEs)
#     else
#         DomainIntegrals.quadrature_d(qs, integrand, element(domain,1), measure, sing)
#     end
# end



doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing = NoSingularity()) =
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, support(dict1, idx1), support(dict2, idx2))

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
            domain1::AbstractInterval, domain2::AbstractInterval)
    integrand = t -> f(t[1],t[2])*unsafe_eval_element(dict1, idx1, t[1])*conj(unsafe_eval_element(dict2,idx2, t[2]))*
        unsafe_weight(measure1, t[1])*unsafe_weight(measure2, t[2])
    DomainIntegrals.integral(integrand, domain1 × domain2, sing)
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicInterval, domain2::PeriodicInterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, d1, d2)
            for d1 in elements(domain1), d2 in elements(domain2))
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicInterval, domain2::UnitInterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, d1, domain2)
            for d1 in elements(domain1))
end


export QuadQBF
"""
QuadQBF uses specialized quadrature routines that incorporate the basis function of the discretization
into their weight function. This means only evaluation of the Green's function are required.
The evaluations are equispaced and can be shared for neighbouring elements in the discretization matrix.
"""
struct QuadQBF{T} <: QuadratureStrategy
    coef    ::  Array{T,1}  #   Coefficients of the two-scale relation
    a       ::  T
    b       ::  T
    x       ::  Array{T,1}  #   1d quadrature points
    w       ::  Array{T,1}  #   1d quadrature weights
end

QuadQBF(args...; options...) = QuadQBF{Float64}(args...; options...)
QuadQBF{T}(n::Int = 1; oversamplingfactor = 2) where {T} =
    QuadQBF{T}(n, oversamplingfactor*(n+1))

function QuadQBF{T}(n::Int, M::Int) where {T}
    coef = bspline_refinable_coefficients(n, T)
    x, w = refinable_quadrature(coef, M)
    QuadQBF{T}(coef, x[1], x[end], x, w)
end

leftpoint(q::QuadQBF) = q.a
rightpoint(q::QuadQBF) = q.b

quad_x(q::QuadQBF) = q.x
quad_w(q::QuadQBF) = q.w



# Map a point x from the interval [a,b] linearly to the interval [c,d]
mapx(x, a, b, c, d) = c + (x-a)/(b-a)*(d-c)


function qbf_quadrature(f, dict, idx, measure, domain, qbf_a, qbf_b, qbf_x, qbf_w)
    a, b = extrema(domain)
    T = typeof(f(a))
    z = zero(T)
    for i = 1:length(qbf_x)
        t = mapx(qbf_x[i], qbf_a, qbf_b, a, b)
        z += qbf_w[i] * f(t) * unsafe_weight(measure, t)
    end
    z * sqrt((b-a) / (qbf_b-qbf_a))
end

function qbf_quadrature2(f, dict1, idx1, measure1, domain1, dict2, idx2, measure2, domain2, qbf_a, qbf_b, qbf_x, qbf_w)
    a1, b1 = extrema(domain1)
    a2, b2 = extrema(domain2)
    T = typeof(f(a1,a2))
    z = zero(T)
    for i = 1:length(qbf_x)
        t1 = mapx(qbf_x[i], qbf_a, qbf_b, a1, b1)
        for j = 1:length(qbf_x)
            t2 = mapx(qbf_x[j], qbf_a, qbf_b, a2, b2)
            z += qbf_w[i] * qbf_w[j] * f(t1,t2) * unsafe_weight(measure1, t1) * unsafe_weight(measure2, t2)
        end
    end
    z * sqrt((b1-a1)*(b2-a2)) / (qbf_b-qbf_a)
end


# We only use QBF quadrature for nonsingular integrals
projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing::NoSingularity, domain::AbstractInterval) =
    qbf_quadrature(f, dict, idx, measure, domain,
        leftpoint(qs), rightpoint(qs), quad_x(qs), quad_w(qs))

# In case of a singularity, we revert to adaptive quadrature
projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing, domain) =
    projectionintegral(QuadAdaptive(), f, dict, idx, measure, sing, domain)

function projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing::NoSingularity, domain::PeriodicInterval)
    if numelements(domain) == 1
        z = projectionintegral(qs, f, dict, idx, measure, sing, element(domain,1))
    else
        A, B = extrema(domain.periodicdomain)
        z = projectionintegral(qs, t -> f(A+mod(t-A,B-A)), dict, idx, measure, sing, domain.subdomain)
    end
    z
end


doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing::NoSingularity,
            domain1::AbstractInterval, domain2::AbstractInterval) =
        qbf_quadrature2(f, dict1, idx1, measure1, domain1, dict2, idx2, measure2, domain2,
            leftpoint(qs), rightpoint(qs), quad_x(qs), quad_w(qs))

function doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing::NoSingularity,
        domain1::PeriodicInterval, domain2::PeriodicInterval)

    A1, B1 = extrema(domain1.periodicdomain)
    A2, B2 = extrema(domain2.periodicdomain)
    f2 = t -> f(A1+mod(t[1]-A1,B1-A1), A2+mod(t[2]-A2,B2-A2))
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, domain1.subdomain, domain2.subdomain)
end



struct QuadGaussLegendre{T} <: QuadratureStrategy
    w   ::  Array{T,1}
    x   ::  Array{T,1}

    function QuadGaussLegendre{T}(n::Int) where {T}
        w, x = gausslegendre(n)
        new(x, w)
    end
end

QuadGaussLegendre(n::Int) = QuadGaussLegendre{Float64}(n::Int)
