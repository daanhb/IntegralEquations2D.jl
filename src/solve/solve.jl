
function solve(Op::BEMOperator, b)
    A = matrix(Op)
    u = A \ b
    Expansion(src(A), u)
end
