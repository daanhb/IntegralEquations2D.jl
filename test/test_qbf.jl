
f_qbftest(t) = 1/(2+cos(2*pi*t))^2

function test_qbf()
    N = 32
    for degree in 1:4
        basis = BSplineTranslatesBasis(N, degree, 0, 1)
        coef_direct = [DomainIntegrals.integral(t -> f_qbftest(t)*BasisFunctions.unsafe_eval_element(basis, i, t), FrameFun.support(basis, i)) for i in 1:length(basis)]
        for oversamplingfactor in 1:5
            A = SimpleIntegralEquations.qbf_operator(basis, degree, oversamplingfactor)
            coef = A * f_qbftest
            @test norm(coef-coef_direct) < 1e-4
        end
    end

    # Do a BigFloat test
    N = 500
    degree = 11
    oversamplingfactor = 6
    basis = BSplineTranslatesBasis(N, degree, big(0.), big(1.))
    A = SimpleIntegralEquations.qbf_operator(basis, degree, oversamplingfactor)
    coef = A * f_qbftest
    # c10 = integral(t -> f_qbftest(t)*BasisFunctions.unsafe_eval_element(basis, 10, t), FrameFun.support(basis, 10))
    c10 = BigFloat("0.004990536452129881961287340961088526093027349985606")
    @test abs(coef[10]-c10) < 1e-45
end
