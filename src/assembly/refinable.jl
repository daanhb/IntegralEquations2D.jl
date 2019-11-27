
"""
Compute the Sweldens quadrature rule for integrals involving a refinable function.

The method is described in:
W. Sweldens and R. Piessens, "Quadrature formulae and asymptotic error expansions
for wavelet approximations of smooth functions", SIAM J. Numer. Anal. 31(4),
pp. 1240-1264, 1994.
"""
function refinable_quadrature(coefficients, M = length(coefficients)-1)
    @assert sum(coefficients) ≈ sqrt(2)

    T = eltype(coefficients)
    L = length(coefficients)
    K = (T(L)-1)/2
    x = collect(range(-K, stop = K, length = M+1))
    w = refinable_quadrature_weights(coefficients, M)
    x, w
end


function refinable_quadrature_weights(coefficients::AbstractVector{T}, M::Int) where T
    # Perform calculations in BigFloat for greater accuracy
    t = range(-BigFloat(1),stop = BigFloat(1), length = M+1)
    A = chebyshev_vandermonde_matrix(t)
    m = chebyshev_moments(BigFloat.(coefficients), M)
    w = A\m
    map(T, w)[:]
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
function chebyshev_moments(coefficients::AbstractVector{T}, M::Int) where {T}
    m = zeros(T,M+1,1)
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
                e = e + T(1)/sqrt(T(2))*coefficients[k]*w[p+1,i+1,k]
            end
            d = d + e*m[i+1]
        end
        m[p+1] = d/(2^p-1)
    end

    m
end

# Calculate the `w' coefficients in the algorithm by Sweldens and Piessens.
# This is based on Appendix A in the above-mentioned reference.
function wcoefficients(T, L, k, M::Int)
    lambda = T(2)*(k-1)/(L-1)-1

    w = zeros(typeof(lambda),M+1,M+1)

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