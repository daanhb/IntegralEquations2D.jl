
using Test

using CompactTranslatesDict,
    DomainIntegrals,
    BasisFunctions,
    FrameFun

using SimpleIntegralEquations

include("test_qbf.jl")

@testset "QBF" begin
    test_qbf()
end
