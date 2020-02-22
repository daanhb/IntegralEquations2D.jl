
function qbf_operator(basis, degree, oversamplingfactor)
    T = domaintype(basis)
    qbf = QuadQBF{T}(degree, oversamplingfactor*(degree+1))
    grid = sampling_grid(basis, oversamplingfactor=oversamplingfactor)
    gb = GridBasis{coefficienttype(basis)}(grid)
    step = oversamplingfactor
    offset = -(length(qbf.x)>>1)
    HorizontalBandedOperator(gb, basis, qbf.w / sqrt(convert(T, length(basis))), step, offset)
end
