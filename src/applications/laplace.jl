"Evaluate the Laplace 2D single layer potential kernel"
function laplace_slp2(x, y)
    T = eltype(x)
    R = norm(x-y)
    if R < eps(T)
         return zero(Complex{T})
    else
        Tpi = convert(T, pi)
        return complex(1/(2*Tpi)*log(R))
    end
end


export Laplace_SLP_2D
"The single layer potential kernel of the Laplace equation."
struct Laplace_SLP_2D{T} <: BoundaryKernel
end

Laplace_SLP_2D() = Laplace_SLP_2D{Float64}()

(kernel::Laplace_SLP_2D)(t, tau, param, x, y) = laplace_slp2(x, y)
(kernel::Laplace_SLP_2D)(x, tau, param, y) = laplace_slp2(x, y)

is_symmetric(::Laplace_SLP_2D) = true

singularity(::Laplace_SLP_2D) = LogSingularDiagonal()

"Evaluate the Laplace 2D double layer potential kernel"
laplace_dlp(x, y, normal_y, z = norm(x-y)) =
    (normal_y' * (x-y)) / z^2/2/pi

"Evaluate the Laplace 2D double layer potential kernel"

function laplace_dlp_kernel(x, y, tau, param)
    R = norm(x-y)
    if R > eps(eltype(x))
        return Complex{Float64}(laplace_dlp(x, y, normal(param, tau), R))
    else
        return Complex{Float64}(0)
    end
end

export Laplace_DLP_2D

"The double layer potential kernel of the Laplace equation in 2D."
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

"Adjoint Double layer potential for Laplace equation"

export Laplace_adjDLP_2D
"The adjoint double layer potential kernel of the Laplace equation in 2D."
struct Laplace_adjDLP_2D{T} <: BoundaryKernel
end

laplace_adjdlp(x, y, normal_x, z = norm(x-y)) =
    (normal_x' * (y-x)) / z^2/2/pi

"Evaluate the adjoint of the Laplace 2D double layer potential kernel"

function laplace_adjdlp_kernel(x, y, t, param)
    z = norm(x-y)
    if abs(z) > 10*eps(eltype(x))
        return Complex{Float64}(laplace_adjdlp(x, y, normal(param, t), z))
    else
        return Complex{Float64}(0)
    end
end

Laplace_adjDLP_2D() = Laplace_adjDLP_2D{Float64}()

BasisFunctions.name(kernel::Laplace_adjDLP_2D) = "Adjoint of 2D Laplace double layer potential kernel"

(kernel::Laplace_adjDLP_2D)(t, tau, param, x, y) =
    laplace_adjdlp_kernel(x, y, t, param)

is_symmetric(::Laplace_adjDLP_2D) = false

# The kernel is continuous but its derivative is singular
singularity(::Laplace_adjDLP_2D) = LogSingularDiagonal()
## Some boundary conditions

export laplace_bcond_harmonic, laplace_param_bcond_harmonic, laplace_parboundary_condition_potential_flow

harmonic_poly(z) = z^2-z+2

boundary_condition_harmonic(x) = real(harmonic_poly(x[1]+im*x[2]))

"Return a field function that evaluates to a plane wave."
laplace_bcond_harmonic() = (x,y) -> boundary_condition_harmonic(SVector(x,y))

"Return a parametric boundary function that evaluates to a plane wave."
laplace_param_bcond_harmonic(param) =
    t -> boundary_condition_harmonic(applymap(param, t))

"Return a parametric boundary function corresponding to the unperturbed flow at r=∞"
laplace_parboundary_condition_potential_flow(param, velocity) =
    t -> velocity'*normal(param, t)
