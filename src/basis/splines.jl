
export qbf_operator

function qbf_operator(basis::Dictionary, degree::Int, oversamplingfactor)
    T = domaintype(basis)
    qbf = QuadQBF{T}(degree, oversamplingfactor*(degree+1))
    _qbf_operator(basis, qbf, oversamplingfactor)
end

function qbf_operator(basis::Dictionary, coef::AbstractVector, oversamplingfactor::Int)
    T = domaintype(basis)
    qbf = QuadQBF{T}(coef; oversamplingfactor=oversamplingfactor)
    _qbf_operator(basis, qbf, oversamplingfactor)
end

function _qbf_operator(basis::Dictionary, qbf::QuadQBF{T}, oversamplingfactor) where {T}
    grid = sampling_grid(basis, oversamplingfactor=oversamplingfactor)
    gb = GridBasis{coefficienttype(basis)}(grid)
    step = oversamplingfactor
    offset = -(length(qbf.x)>>1)
    HorizontalBandedOperator(gb, basis, qbf.w / sqrt(convert(T, length(basis))), step, offset)
end
