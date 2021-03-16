
"""
Compute the Sweldens quadrature rule for integrals involving a refinable function.

The method is described in:
W. Sweldens and R. Piessens, "Quadrature formulae and asymptotic error expansions
for wavelet approximations of smooth functions", SIAM J. Numer. Anal. 31(4),
pp. 1240-1264, 1994.

We return a quadrature rule that with points spaced 1/M apart, symmetric with
respect to the origin, and including the origin itself.
"""
function refinable_quadrature(coefficients, M = length(coefficients)-1)
    @assert sum(coefficients) ≈ sqrt(2)
    L = length(coefficients)
    # Perform calculations in BigFloat for greater accuracy
    if isodd(L) || (iseven(L) && iseven(M))
        t = range(-BigFloat(1.),stop = BigFloat(1.), length = M+1)
    else
        # We have to shift the grid in order to include the center point
        t = range(-BigFloat(1.),stop = BigFloat(1.), length = 2M+1)[2:2:end-1]
    end
    A = chebyshev_vandermonde_matrix(t)
    m = chebyshev_moments(BigFloat.(coefficients), length(t)-1)
    w = A\m
    (L-1) * t/2, map(eltype(coefficients), w)[:]
end


function chebyshev_vandermonde_matrix(t::AbstractVector{T}) where {T}
    M = length(t)-1
    A = zeros(T, M+1, M+1)
    for i in 0:M
        A[i+1,:] = chebyshev_eval.(i, t)
    end
    A
end

"Evaluate the chebyshev polynomial of degree p in x."
chebyshev_eval(p, x) = real(cos(p*acos(x)))


"""
Calculate the modified moments of the scaling function defined by coefficients.
"""
function chebyshev_moments(coefficients::AbstractVector, M::Int)
    T = eltype(coefficients)
    m = zeros(T,M+1)
    m[1] = 1
    L = length(coefficients)

    w = zeros(T,M+1,M+1,L)
    for k in 1:L
        w[:,:,k] = wcoefficients(T, L, k, M)
    end

    for p in 1:M
        d = zero(T)
        for i in 0:p-1
            e = zero(T)
            for k in 1:L
                # % CAVE: our coefficients differ by a factor 1/sqrt(2) from Sweldens!
                e = e + 1/sqrt(T(2))*coefficients[k]*w[p+1,i+1,k]
            end
            d = d + e*m[i+1]
        end
        m[p+1] = d/(2^p-1)
    end

    m
end

# Calculate the `w' coefficients in the algorithm by Sweldens and Piessens.
# This is based on Appendix A in the above-mentioned reference.
function wcoefficients(::Type{T}, L, k, M::Int) where {T}
    lambda = T(2)*(k-1)/(T(L)-1)-1

    w = zeros(T,M+1,M+1)

    w[1,1] = 1
    w[2,1] = lambda
    w[2,2] = 1
    if M > 1
        w[3,1] = 2*lambda^2-3
        w[3,2] = 4*lambda
        w[3,3] = 1
    end

    # p+1 is the degree of the next polynomial
    for p in 2:M-1
        w[p+2,1] = w[p+1,2] + 2*lambda*w[p+1,1] - 4*w[p,1]
        w[p+2,2] = 2*w[p+1,1]+w[p+1,3]+2*lambda*w[p+1,2]-4*w[p,2]
        for i in 2:p-1
            w[p+2,i+1] = w[p+1,i]+w[p+1,i+2]+2*lambda*w[p+1,i+1]-4*w[p,i+1]
        end
        w[p+2,p+1] = w[p+1,p] + 2*lambda*w[p+1,p+1]
        w[p+2,p+2] = w[p+1,p+1]
    end

    w
end


bspline_refinable_coefficients(degree::Int, ::Type{T} = Float64) where {T} =
    [binomial(degree+1,k) for k in 0:degree+1] / 2^(degree+one(T)/2)



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

QuadQBF{T}(n::Int, M::Int) where {T} = QuadQBF{T}(bspline_refinable_coefficients(n, T), M)

function QuadQBF{T}(coef::AbstractVector, M::Int) where {T}
    x, w = refinable_quadrature(coef, M)
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
