
using Test

using BasisFunctions,
    BasisTranslates, BasisTranslates.BSplines,
    DomainIntegrals,
    DomainSets,
    FrameFun,
    StaticArrays

using IntegralEquations2D

include("test_qbf.jl")
include("test_helmholtz.jl")
include("test_laplace.jl")

# @testset "QBF" begin
#     test_qbf()
# end

@testset "Helmholtz BEM" begin
    test_helmholtz()
end
@testset "Laplace BEM" begin
    test_laplace()
end
