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



## Some boundary conditions

export laplace_bcond_harmonic, laplace_param_bcond_harmonic

harmonic_poly(z) = z^2-z+2

boundary_condition_harmonic(x) = real(harmonic_poly(x[1]+im*x[2]))

"Return a field function that evaluates to a plane wave."
laplace_bcond_harmonic() = (x,y) -> boundary_condition_harmonic(SVector(x,y))

"Return a parametric boundary function that evaluates to a plane wave."
laplace_param_bcond_harmonic(param) =
    t -> boundary_condition_harmonic(applymap(param, t))
