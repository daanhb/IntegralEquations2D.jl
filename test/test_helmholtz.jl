IE = SimpleIntegralEquations

function test_helmholtz()
    wavenumber = 5.0
    N = 128
    M = N
    splinedegree = 1
    obstacle = IE.Kite(2)
    direction = SVector(1.0, 0.0)
    amplitude = 1.0
    param = parameterization(obstacle)
    paramdomain = domain(param)
    SingleLayerPotential = Helmholtz_SLP_2D(wavenumber)
    BIO = BoundaryIntegralOperator(SingleLayerPotential, obstacle, param)
    basis = complex(BSplineTranslatesBasis(N, splinedegree, leftendpoint(paramdomain), rightendpoint(paramdomain)))
    basis_obstacle = BasisFunctions.ParamDict(basis, param, obstacle)
    coll_points = PeriodicEquispacedGrid(M, paramdomain)
    coll_points_obstacle = mapped_grid(coll_points, param)
    sampling_col = GridSampling(coll_points, Complex{Float64})
    sampling_col_obstacle = GridSampling(coll_points_obstacle, Complex{Float64})
    sampling_gal = ProjectionSampling(complex(basis), IE.measure(BIO))
    bcond = make_parboundary_condition_planewave(param, wavenumber, direction, amplitude)
    bcond_field = make_boundary_condition_planewave(wavenumber, direction, amplitude)
    BEM_col = (sampling_col * BIO) * basis
    BEM_gal = (sampling_gal * BIO) * basis
    quad_qbf = QuadQBF(splinedegree)
    quad_gk = QuadAdaptive()

    z1 = compute_BEM_entry(BEM_col, 2, 2, quad_qbf)
    @test abs(z1 - (0.2296280626495274 + 0.11152548885226357im)) < 1e-3

    IE.assemble!(BEM_col, quad_qbf)
    IE.assemble!(BEM_gal, quad_qbf)
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
    @test abs(z_col-z_exact) / abs(z_exact) < 1e-4

    z_gal = eval_field(BEM_gal, density_gal, point)
    @test abs(z_gal - z_exact) / abs(z_exact) < 1e-5
end
