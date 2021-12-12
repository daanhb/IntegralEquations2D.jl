"Evaluate the Laplace 2D single layer potential kernel"
function laplace_slp(x, y)
    T = eltype(x)
    if norm(x-y) < eps(T)
         return zero(Complex{T})
    else
        return complex(log(norm(x-y)))
    end
end

export Laplace_SLP_2D
"The single layer potential kernel of the Laplace equation."
struct Laplace_SLP_2D{T} <: BoundaryKernel
end

Laplace_SLP_2D() = Laplace_SLP_2D{Float64}()

(kernel::Laplace_SLP_2D)(t, tau, param, x, y) = laplace_slp(x, y)
(kernel::Laplace_SLP_2D)(x, tau, param, y) = laplace_slp(x, y)

is_symmetric(::Laplace_SLP_2D) = true

singularity(::Laplace_SLP_2D) = LogSingularDiagonal()


"Evaluate the Laplace 2D double layer potential kernel"
hh_dlp(x, y, wavenumber, normal_y, z = norm(x-y)) =
    (normal_y' * (x-y)) / z^2


"Evaluate the Laplace 2D double layer potential kernel"
function laplace_dlp_kernel(x, y, tau, wavenumber, param)
    z = norm(x-y)
    if abs(z) > eps(eltype(x))
        complex(laplace_dlp(x, y, wavenumber, normal(param, tau), z))
    else
        return zero(eltype(x))
    end
end

export Laplace_DLP_2D
"The double layer potential kernel of the Helmholtz equation in 2D."
struct Laplace_DLP_2D{T} <: BoundaryKernel
end

Laplace_DLP_2D() = Laplace_DLP_2D{Float64}()

BasisFunctions.name(kernel::Laplace_DLP_2D) = "2D Laplace double layer potential kernel"

(kernel::Laplace_DLP_2D)(t, tau, param, x, y) =
    laplace_dlp_kernel(x, y, tau, param)
(kernel::Laplace_DLP_2D)(x, tau, param, y) =
    laplace_dlp_kernel(x, y, tau, param)

is_symmetric(::Laplace_DLP_2D) = false

# The kernel is continuous but its derivative is singular
singularity(::Laplace_DLP_2D) = LogSingularDiagonal()


## Some boundary conditions

export laplace_bcond_harmonic, laplace_param_bcond_harmonic

harmonic_poly(z) = z^2-z+2

boundary_condition_harmonic(x) = real(harmonic_poly(x[1]+im*x[2]))

"Return a field function that evaluates to a plane wave."
laplace_bcond_harmonic() = (x,y) -> boundary_condition_harmonic(SVector(x,y))

"Return a parametric boundary function that evaluates to a plane wave."
laplace_param_bcond_harmonic(param) =
    t -> boundary_condition_harmonic(applymap(param, t))
