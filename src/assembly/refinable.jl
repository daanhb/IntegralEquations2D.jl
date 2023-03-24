
################################
# Definition of QBF quadrature
################################

export QuadQBF
"""
QuadQBF uses specialized quadrature routines that incorporate the basis function of the discretization
into their weight function. This means only evaluation of the Green's function are required.
The evaluations are equispaced and may be shared for neighbouring elements in the discretization matrix.
"""
struct QuadQBF{T} <: DomainIntegrals.IntervalRule{T}
    coef    ::  Array{T,1}  #   Coefficients of the two-scale relation
    a       ::  T
    b       ::  T
    x       ::  Array{T,1}  #   1d quadrature points
    w       ::  Array{T,1}  #   1d quadrature weights
end

QuadQBF(args...; options...) = QuadQBF{Float64}(args...; options...)
QuadQBF(coef::AbstractArray{T}, args...; options...) where {T} = QuadQBF{T}(coef, args...; options...)
QuadQBF{T}(n::Int = 1; oversamplingfactor = 2) where {T} =
    QuadQBF{T}(n, oversamplingfactor*(n+1))
QuadQBF{T}(coef::AbstractVector; oversamplingfactor = 2) where {T} =
    QuadQBF{T}(coef, oversamplingfactor*length(coef)-1)

QuadQBF{T}(n::Int, M::Int) where {T} =
    QuadQBF{T}(RefinableFunctions.bspline_refinable_coefficients(n, T), M)

function QuadQBF{T}(coef::AbstractVector, M::Int) where {T}
    x, w = RefinableFunctions.refinable_quadrature(coef, M)
    # shift the quadrature points to make them symmetric
    x = x .- (length(coef)-1)/2
    QuadQBF{T}(coef, -(length(coef)-one(T))/2, (length(coef)-one(T))/2, x, w)
end

leftpoint(q::QuadQBF) = q.a
rightpoint(q::QuadQBF) = q.b

points(q::QuadQBF) = q.x
weights(q::QuadQBF) = q.w

DomainIntegrals.domain(q::QuadQBF) = q.a..q.b

splineorder(q::QuadQBF) = length(q.coef)-2



# Map a point x from the interval [a,b] linearly to the interval [c,d]
mapx(x, a, b, c, d) = c + (x-a)/(b-a)*(d-c)

function qbf_quadrature(f, dict, idx, measure, domain, qbf_a, qbf_b, qbf_x, qbf_w)
    a, b = extrema(domain)
    T = typeof(f(a))
    z = zero(T)
    for i = 1:length(qbf_x)
        t = mapx(qbf_x[i], qbf_a, qbf_b, a, b)
        z += qbf_w[i] * f(t) * unsafe_weightfun(measure, t)
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
            z += qbf_w[i] * qbf_w[j] * f(t2,t1) * unsafe_weightfun(measure1, t1) * unsafe_weightfun(measure2, t2)
        end
    end
    z * sqrt((b1-a1)*(b2-a2)) / (qbf_b-qbf_a)
end
