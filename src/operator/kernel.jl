
"A kernel is a function of two variables."
abstract type Kernel end

"A boundary kernel is associated with a parameterized boundary."
abstract type BoundaryKernel <: Kernel end

# A concrete kernel k should implement the routine k(x,y,param).
# The boundary kernel may not know its own parameterization, but
# the parameterization is passed as an extra argument.

# By default, we add a parameterization by asking the kernel
(kernel::BoundaryKernel)(x, y) = kernel(x, y, parameterization(kernel))
# If two numbers are given, we assume they are in the parameter domain
(kernel::BoundaryKernel)(t::Number, tau::Number, param) =
    kernel(t, tau, param, applymap(param, t), applymap(param, tau))
(kernel::BoundaryKernel)(x::SVector, tau::Number, param) =
    kernel(x, tau, param, applymap(param, tau))
