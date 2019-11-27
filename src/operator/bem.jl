
"The abstract supertype of all BEM matrix operators."
abstract type BEMOperator{T} <: DictionaryOperator{T} end

"A dense BEM operator stores the dense matrix."
mutable struct DenseBEMOperator{S,T} <: BEMOperator{T}
    src         ::  Dictionary
    sampling    ::  SamplingOperator
    intop       ::  IntegralOperator{S,T}
    quad        ::  QuadratureStrategy
    A           ::  Array{T,2}
    isassembled ::  Bool
end

initvalue(::Type{T}) where {T <: AbstractFloat} = floatmax(T)
initvalue(::Type{Complex{T}}) where {T <: AbstractFloat} = initvalue(real(T))

function DenseBEMOperator(dict::Dictionary, sampling::SamplingOperator,
        intop::IntegralOperator{S,T}, quad = QuadAdaptive()) where {S,T}
    A = fill(initvalue(T), size(sampling)[1],length(dict))
    DenseBEMOperator{S,T}(dict, sampling, intop, quad, A, false)
end

integraloperator(A::DenseBEMOperator) = A.intop

BasisFunctions.dest(op::DenseBEMOperator) = dest(op.sampling)

isassembled(A::DenseBEMOperator{S,T}, i, j) where {S,T} = !(A.A[i,j] == initvalue(T))
isassembled(A::DenseBEMOperator) = A.isassembled == true

Base.:*(op::DiscretizedIntOp, dict::Dictionary) =
    DenseBEMOperator(dict, op.sampling, op.intop)

assemble!(A::DenseBEMOperator, quad = A.quad; verbose = false) =
    _assemble!(A, quad, A.src, A.sampling, A.intop; verbose = verbose)

function _assemble!(A::DenseBEMOperator, quad, dict, sampling, intop; verbose = false)
    for i in 1:size(sampling)[1]
        verbose && println("Assembly: row $i")
        for j in 1:length(dict)
            A.A[i,j] = compute_BEM_entry(A, i, j, quad)
        end
    end
    A.isassembled = true
    A.A
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
