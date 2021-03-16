
besselh0(z) = besselj0(z) + im*bessely0(z)
besselh1(z) = besselj1(z) + im*bessely1(z)

"Evaluate the Helmholtz 2D single layer potential kernel"
function hh_slp(x, y, wavenumber)
    T = eltype(x)
    if norm(x-y) < eps(T)
        return zero(Complex{T})
    else
        return im * besselh0(wavenumber*norm(x-y)) / 4
    end
end

"Evaluate the Helmholtz 2D double layer potential kernel"
hh_dlp(x, y, wavenumber, normal_y, z = norm(x-y)) =
    -(normal_y' * (x-y)) * im * wavenumber/4 * besselh1(wavenumber*z) / z

"Diagonal limiting value of the Helmholtz 2D double layer potential kernel"
function hh_dlp_limit(x, y, wavenumber, normal_y, grad_y, grad_d_y)
    -(grad_d_y[1] * grad_y[2] - grad_y[1] * grad_d_y[2]) / (4*norm(grad_y)^3*pi) +0im
end

"Evaluate the Helmholtz 2D double layer potential kernel"
function hh_dlp_kernel(x, y, tau, wavenumber, param)
    z = norm(x-y)
    if abs(z) > 100*eps(eltype(x))
        hh_dlp(x, y, wavenumber, normal(param, tau), z)
    else
        # Evaluate the limiting value
        hh_dlp_limit(x, y, wavenumber, normal(param, tau),
            gradient(param, tau), gradient_derivative(param, tau))
    end
end

hh_slpt(x, y, wavenumber, normal_x) = error("to do")
hh_hypersingular(x, y, wavenumber, normal_x, normal_y) = error("to do")



"Supertype of all Helmholtz kernels."
abstract type HelmholtzKernel <: BoundaryKernel end

"Wavenumber of the Helmholtz equation"
wavenumber(op::HelmholtzKernel) = op.wavenumber

export Helmholtz_SLP_2D
"The single layer potential kernel of the Helmholtz equation."
struct Helmholtz_SLP_2D{T} <: HelmholtzKernel
    wavenumber  ::  T
end

Helmholtz_SLP_2D(wavenumber::Integer) = Helmholtz_SLP_2D(float(wavenumber))

BasisFunctions.name(kernel::Helmholtz_SLP_2D) = "2D Helmholtz single layer potential kernel"

# Two arguments are given in the parameter domain
(kernel::Helmholtz_SLP_2D)(t::Number, tau::Number, param) =
    kernel(t, tau, param, applymap(param, t), applymap(param, tau))
(kernel::Helmholtz_SLP_2D)(t, tau, param, x, y) = hh_slp(x, y, wavenumber(kernel))

# Or: the first argument is a field point
(kernel::Helmholtz_SLP_2D)(x::SVector{2}, tau::Number, param) =
    kernel(x, tau, param, applymap(param, tau))
(kernel::Helmholtz_SLP_2D)(x, tau, param, y) = hh_slp(x, y, wavenumber(kernel))

is_symmetric(::Helmholtz_SLP_2D) = true

singularity(::Helmholtz_SLP_2D) = LogDiagonallySingular()


export Helmholtz_DLP_2D
"The double layer potential kernel of the Helmholtz equation in 2D."
struct Helmholtz_DLP_2D{T} <: HelmholtzKernel
    wavenumber  ::  T
end

Helmholtz_DLP_2D(wavenumber::Integer) = Helmholtz_DLP_2D(float(wavenumber))

BasisFunctions.name(kernel::Helmholtz_DLP_2D) = "2D Helmholtz double layer potential kernel"

(kernel::Helmholtz_DLP_2D)(t::Number, tau::Number, param) =
    kernel(t, tau, param, applymap(param, t), applymap(param, tau))
(kernel::Helmholtz_DLP_2D)(t, tau, param, x, y) =
    hh_dlp_kernel(x, y, tau, wavenumber(kernel), param)

(kernel::Helmholtz_DLP_2D)(x::SVector{2}, tau::Number, param) =
    kernel(x, tau, param, applymap(param, tau))
(kernel::Helmholtz_DLP_2D)(x, tau, param, y) =
    hh_dlp_kernel(x, y, tau, wavenumber(kernel), param)

is_symmetric(::Helmholtz_DLP_2D) = false

# The kernel is continuous but its derivative is singular
singularity(::Helmholtz_DLP_2D) = LogDiagonallySingular()
