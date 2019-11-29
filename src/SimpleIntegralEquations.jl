module SimpleIntegralEquations

using QuadGK, HCubature, FastGaussQuadrature,
    SpecialFunctions, StaticArrays,
    IntervalSets, DomainSets,
    DomainIntegrals,
    BasisFunctions, GridArrays,
    CardinalBSplines, CompactTranslatesDict,
    LinearAlgebra,
    RecipesBase


import DomainSets:
    Domain,
    EuclideanDomain,
    AbstractMap,
    applymap,
    jacobian,
    gradient,
    parameterization,
    domain

import BasisFunctions:
    AbstractOperator,
    GridSampling,
    ProjectionSampling,
    SamplingOperator,
    measure,
    integral,
    period,
    matrix

import DomainIntegrals:
    Singularity,
    NoSingularity,
    CurveSingularity,
    DiagonalSingularity,
    LogPointSingularity

export
    measure

const Vec2{T} = SVector{2,T}
const Vec3{T} = SVector{3,T}

include("kernels/kernel.jl")
include("kernels/helmholtz.jl")

include("domains/param.jl")
include("domains/smooth.jl")
include("domains/recipes.jl")

include("assembly/refinable.jl")
include("assembly/quadrature.jl")
include("assembly/sampling.jl")

include("operator/operator.jl")
include("operator/bem.jl")
include("operator/bie.jl")
include("operator/bndcondition.jl")

include("assembly/assembly.jl")
include("assembly/field.jl")

include("solve/solve.jl")

end # module
