
f_qbftest(t) = 1/(2+cos(2*pi*t))^2

function test_qbf()
    N = 32
    for degree in 1:4
        basis = 2^(log2(N)/2)*PeriodicBSplines(N; degree)
        coef_direct = [DomainIntegrals.integral(t -> f_qbftest(t)*BasisFunctions.unsafe_eval_element(basis, i, t), FrameFun.support(basis, i)) for i in 1:length(basis)]
        for oversamplingfactor in 1:5
            A = IntegralEquations2D.qbf_operator(basis, degree, oversamplingfactor)
            coef = A * f_qbftest
            @test norm(coef-coef_direct) < 1e-4
        end
    end

    # Do a BigFloat test
    N = 500
    degree = 11
    oversamplingfactor = 6
    basis = 2^(log2(big(N))/2)*PeriodicBSplines{BigFloat}(N; degree=degree)
    A = IntegralEquations2D.qbf_operator(basis, degree, oversamplingfactor)
    coef = A * f_qbftest
    # c10 = integral(t -> f_qbftest(t)*BasisFunctions.unsafe_eval_element(basis, 10, t), FrameFun.support(basis, 10))
    c10 = BigFloat("0.004990536452129881961287340961088526093027349985606")
    @test abs(coef[10]-c10) < 1e-45
end
