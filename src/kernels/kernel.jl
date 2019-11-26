
"""
A kernel is a function `K(x,y)` of two variables, where `x` and `y` have type
`S` and the result has type `T`.
"""
abstract type Kernel{S,T} end

BasisFunctions.name(kernel::Kernel) = "Generic kernel function"

isreal(k::Kernel{S,T}) where {S,T} = false
isreal(k::Kernel{S,T}) where {S,T <: Real} = true

# Logarithmic singularity of 2D Green's functions
struct LogDiagonalSingularity <: CurveSingularity end
# The typical 1/R singularity
struct RadialDiagonalSingularity <: CurveSingularity end
# Generic, unknown singularity
struct UnknownDiagonalSingularity <: CurveSingularity end

# By default, we know nothing
singularity(kernel::Kernel) = UnknownSingularity()
