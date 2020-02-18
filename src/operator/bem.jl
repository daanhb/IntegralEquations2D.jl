
"The abstract supertype of all BEM matrix operators."
abstract type BEMOperator{T} <: DictionaryOperator{T} end

kernel(A::BEMOperator) = kernel(integraloperator(A))
kernelfunction(A::BEMOperator) = kernelfunction(integraloperator(A))
measure(A::BEMOperator) = measure(integraloperator(A))
singularity(A::BEMOperator) = singularity(integraloperator(A))


"A dense BEM operator stores the dense matrix."
mutable struct DenseBEMOperator{T} <: BEMOperator{T}
    src         ::  Dictionary
    sampling    ::  SamplingOperator
    intop       ::  IntegralOperator
    quad        ::  QuadratureStrategy
    A           ::  Array{T,2}
    isassembled ::  Bool
end

initvalue(::Type{T}) where {T <: AbstractFloat} = floatmax(T)
initvalue(::Type{Complex{T}}) where {T <: AbstractFloat} = initvalue(real(T))

function DenseBEMOperator(dict::Dictionary, sampling::SamplingOperator,
        intop::IntegralOperator, quad = QuadAdaptive())
    T = promote_type(coefficienttype(dict), codomaintype(dest(sampling)))
    A = fill(initvalue(T), size(sampling)[1],length(dict))
    DenseBEMOperator{T}(dict, sampling, intop, quad, A, false)
end

samplingoperator(A::DenseBEMOperator) = A.sampling
integraloperator(A::DenseBEMOperator) = A.intop


BasisFunctions.dest(op::DenseBEMOperator) = dest(op.sampling)

isassembled(A::DenseBEMOperator{T}, i, j) where {T} =
    !(A.A[i,j] == initvalue(T))
isassembled(A::DenseBEMOperator) = A.isassembled == true

Base.:*(op::SampledIntegralOperator, dict::Dictionary) =
    DenseBEMOperator(dict, samplingoperator(op), integraloperator(op))

export assemble!
"Assemble the dense BEM matrix."
function assemble!(op::DenseBEMOperator, quad = op.quad, A = op.A; verbose = false)
    for i in 1:size(A, 1)
        verbose && println("Assembly: row $i")
        for j in 1:size(A,2)
            A[i,j] = compute_BEM_entry(op, i, j, quad)
        end
    end
    op.isassembled = true
    A
end

BasisFunctions.matrix(A::DenseBEMOperator) =
    isassembled(A) ? A.A : assemble!(A)

# We cache the computed entry, because why not, the memory is allocated anyway
Base.getindex(A::DenseBEMOperator, i::Int, j::Int) =
    isassembled(A, i, j) ? A.A[i,j] : A.A[i,j] = compute_BEM_entry(A, i, j)

function Base.getindex(A::DenseBEMOperator, i, j)
    if !isassembled(A)
        compute_BEM_range!(A.A, A, i, j)
    end
    A.A[i,j]
end
