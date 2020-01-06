
"""
Supertype of all Helmholtz kernels.

They all support a wavenumber. They map an argument in Euclidean space to
a complex number.
"""
abstract type HelmholtzKernel{N,T} <: Kernel{SVector{N,T},Complex{T}} end

wavenumber(k::HelmholtzKernel) = k.wavenumber


"""
The single layer potential kernel of the Helmholtz equation in 2D.
"""
struct Helmholtz_SLP_2D{S,T} <: HelmholtzKernel{2,T}
    wavenumber  ::  S
end

BasisFunctions.name(kernel::Helmholtz_SLP_2D) = "2D Helmholtz single layer potential kernel"

Helmholtz_SLP_2D(wavenumber::Integer) = Helmholtz_SLP_2D(float(wavenumber))

Helmholtz_SLP_2D(wavenumber::T) where {T<:Real} =
    Helmholtz_SLP_2D{T,T}(wavenumber)
Helmholtz_SLP_2D(wavenumber::T) where {T<:Complex} =
    Helmholtz_SLP_2D{T,real(T)}(wavenumber)

is_symmetric(k::Helmholtz_SLP_2D) = true

singularity(k::Helmholtz_SLP_2D) = LogDiagonalSingularity()


# """
# Evaluate the Hankel function of the first kind and of order nu from its
# asymptotic form for large arguments.
#
# The implementation is based on Abramowitz and Stegun, (9.2.7).
# """
# function besselh_asymptotic(nu, z::Real)
#     mu = 4*nu^2
#     P = 1 - (mu-1)*(mu-9)*(8*z)^(-2)/factorial(2) +
#         (mu-1)*(mu-9)*(mu-25)*(mu-49)*(8*z)^(-4)/factorial(4)
#     Q = (mu-1)*(8*z)^(-1) - (mu-1)*(mu-9)*(mu-25)*(8*z)^(-3)/factorial(3)
#     chi = z - (1/2*nu+1/4)*pi
#     sqrt(2*z^(-1)/pi) * (P + im*Q) * exp(im*chi)
# end
#
# function besselh_asymptotic(nu, z::Complex)
#     mu = 4*nu^2
#     P = 1 - (mu-1)*(mu-9)*(8*z)^(-2)/factorial(2) +
#         (mu-1)*(mu-9)*(mu-25)*(mu-49)*(8*z)^(-4)/factorial(4)
#     Q = (mu-1)*(8*z)^(-1) - (mu-1)*(mu-9)*(mu-25)*(8*z)^(-3)/factorial(3)
#
#     chi = real(z) - (1/2*nu+1/4)*pi
#     h = sqrt(2*z^(-1)/pi) * (P + im*Q) * exp(im*chi)
#     h * exp(-imag(z)) # In case of overflow or underflow, this is the likely culprit
# end


"Assume that `K(x,y) = phi(z)` with `z=|x-y|`."
function eval_phi(kernel::Helmholtz_SLP_2D, z)
    kz = wavenumber(kernel) * z
    (im*besselj(0, kz) - bessely(0, kz)) / 4 # AMOS' zbesh (besselh) is about 1.5 times slower
end
# eval_phi(kernel::Helmholtz_SLP_2D, z) = im * besselh(0, 1, wavenumber(kernel)*z) / 4

eval_kernel(kernel::Helmholtz_SLP_2D, x, y) = eval_phi(kernel, norm(x-y))


kernel_singular_reg(kernel::Helmholtz_SLP_2D{T}, x::Vec2{T}, y::Vec2{T}) where {T} =
    -1/4 * 2/T(pi) * besselj(0, kernel.wavenumber * norm(x-y))


"""
The double layer potential kernel of the Helmholtz equation in 2D.
"""
struct Helmholtz_DLP_2D{U,S,T} <: HelmholtzKernel{S,T}
    wavenumber  ::  U
end

BasisFunctions.name(kernel::Helmholtz_DLP_2D) = "2D Helmholtz double layer potential kernel"

# kernel_singular_reg(kernel::Helmholtz_SLP_2D{T}, param::Boundary{T}, t::T, tau::T) where {T} =
#     kernel_singular_reg(kernel, param(t), param(tau))
#
#
# function kernel_regular(kernel::Helmholtz_SLP_2D{T}, param::Boundary{T}, t::T, tau::T, N) where {T}
#     x = param(t)
#     y = param(tau)
#     res = norm(x-y)
#
#     if res > 10eps(T)
#         kern = kernel(x, y)
#         approx = log(abs(t-tau)) * kernel_singular_reg(kernel, param, t, tau)
#         z1 = kern - approx
#     else
#         gamma = 0.5772156649015328606*one(T)
#         grad = gradient(param, t)
#         rt = kernel.wavenumber * norm(grad)
#         z1 =  im/4 - 1/4*2/T(pi)*(log(rt/2)+gamma)
#     end
#
#     z1 + 1/4*2/pi*log(N/period(param))*besselj(0, kernel.wavenumber * res)
# end
