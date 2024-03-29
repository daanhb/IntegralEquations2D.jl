
function test_helmholtz(T = Float64, N = 128)
    wavenumber = T(5)
    M = 2N
    splinedegree = 1
    obstacle = Kite{T}(2)
    direction = SVector(one(T), zero(T))
    amplitude = one(T)
    param = parameterization(obstacle)
    paramdomain = domain(param)
    SingleLayerPotential = Helmholtz_SLP_2D(wavenumber)
    BIO = BoundaryIntegralOperator(SingleLayerPotential, obstacle, param)
    if paramdomain == Interval{:closed,:open,T}(0,1)
        basis = complex(2^(log2(T(N))/2)*PeriodicBSplines{T}(N, degree=splinedegree))
    else
        a, b = extrema(paramdomain)
        m = mapto(a..b, 0..1)
        splines = 2^(log2(T(N))/2)*PeriodicBSplines{T}(N, degree=splinedegree)
        basis = complex(MappedDict(splines, m))
    end
    basis_obstacle = BasisFunctions.ParamDict(basis, param, obstacle)
    coll_points = PeriodicEquispacedGrid(M, paramdomain)
    coll_points_obstacle = map_grid(param, coll_points)
    sampling_col = GridSampling(coll_points, Complex{T})
    sampling_col_obstacle = GridSampling(coll_points_obstacle, Complex{T})
    sampling_gal = ProjectionSampling(basis, measure(BIO))
    bcond = make_parboundary_condition_planewave(param, wavenumber, direction, amplitude)
    bcond_field = make_boundary_condition_planewave(wavenumber, direction, amplitude)
    BEM_col = (sampling_col * BIO) * basis
    BEM_gal = (sampling_gal * BIO) * basis
    quad_qbf = QuadQBF(splinedegree; oversamplingfactor=2)
    quad_adaptive = QuadAdaptive()
    bemquad_qbf_graded = BEMQuadQBF_Graded(quad_qbf)
    bemquad_qbf_adaptive = BEMQuadQBF_Adaptive(quad_qbf, quad_adaptive)

    z1 = compute_BEM_entry(BEM_col, 2, 2, bemquad_qbf_graded)
    @test abs(z1 - (0.19220429036977982 + 0.11128991068138182im)) < 1e-3

    assemble!(BEM_col, bemquad_qbf_graded)
    assemble!(BEM_gal, bemquad_qbf_adaptive)
    A_col = copy(matrix(BEM_col))
    b_col = sampling_col * bcond
    coef_col = A_col \ b_col
    density_col = Expansion(basis, coef_col)
    A_gal = copy(matrix(BEM_gal))
    b_gal = sampling_gal * bcond
    coef_gal = A_gal \ b_gal
    density_gal = Expansion(basis, coef_gal)

    point = SVector(0.05, -0.2)
    z_exact = bcond_field(point...)

    z_col = eval_field(BEM_col, density_col, point)
    @test abs(z_col-z_exact) / abs(z_exact) < 1e-3

    z_gal = eval_field(BEM_gal, density_gal, point)
    @test abs(z_gal - z_exact) / abs(z_exact) < 1e-5
end
