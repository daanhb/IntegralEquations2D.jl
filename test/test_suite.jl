
using Test

using BasisFunctions,
    CompactTranslatesDict,
    DomainIntegrals,
    DomainSets,
    FrameFun,
    StaticArrays

using SimpleIntegralEquations

include("test_qbf.jl")
include("test_helmholtz.jl")

# @testset "QBF" begin
#     test_qbf()
# end

@testset "Helmholtz BEM" begin
    test_helmholtz()
end
