module IntegralEquations2D

using QuadGK, HCubature, FastGaussQuadrature,
    SpecialFunctions, StaticArrays,
    IntervalSets, DomainSets,
    DomainIntegrals,
    BasisFunctions, FrameFun, GridArrays,
    BasisTranslates, BasisTranslates.BSplines,
    LinearAlgebra,
    RecipesBase

import Base:
    Threads

import DomainSets:
    Domain,
    EuclideanDomain,
    applymap,
    jacobian,
    parameterization

import BasisFunctions:
    AbstractOperator,
    GridSampling,
    ProjectionSampling,
    SamplingOperator,
    measure,
    integral,
    period,
    matrix,
    domain

import DomainIntegrals:
    Singularity,
    NoSingularity,
    CurveSingularity,
    SingularDiagonal,
    LogSingPoint

export measure

const Vec2{T} = SVector{2,T}
const Vec3{T} = SVector{3,T}


include("util/common.jl")

include("domains/param.jl")
include("domains/polygons.jl")
include("domains/smooth.jl")
include("domains/recipes.jl")

include("assembly/refinable.jl")
include("assembly/quadrature.jl")
include("assembly/sampling.jl")

include("basis/splines.jl")

include("operator/kernel.jl")
include("operator/operator.jl")
include("operator/bie.jl")
include("operator/bem.jl")

include("assembly/assembly.jl")
include("assembly/field.jl")

include("solve/solve.jl")

include("applications/helmholtz.jl")
include("applications/laplace.jl")

end # module
