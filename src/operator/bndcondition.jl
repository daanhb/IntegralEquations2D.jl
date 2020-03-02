
export make_boundary_condition_planewave,
    make_parboundary_condition_planewave,
    make_boundary_condition_pointsource,
    make_parboundary_condition_pointsource

boundary_condition_wave(x, wavenumber, direction, amplitude) =
    amplitude * exp(im*wavenumber*dot(direction, x))

"Return a field function that evaluates to a plane wave."
make_boundary_condition_planewave(wavenumber, direction, amplitude = 1) =
    (x,y) -> boundary_condition_wave(SVector(x,y), wavenumber, direction, amplitude)

"Return a parametric boundary function that evaluates to a plane wave."
make_parboundary_condition_planewave(param, wavenumber, direction, amplitude = 1) =
    t -> boundary_condition_wave(applymap(param, t), wavenumber, direction, amplitude)

boundary_condition_pointsource(x, center, wavenumber, amplitude) =
    amplitude * besselh(0, 1, wavenumber*norm(x-center))

"Return a field function that evaluates to a point source."
make_boundary_condition_pointsource(center, wavenumber, amplitude = 1) =
    (x,y) -> boundary_condition_pointsource(SVector(x,y), center, wavenumber, amplitude)

"Return a parametric boundary function that evaluates to a point source."
make_parboundary_condition_pointsource(param, center, wavenumber, amplitude = 1) =
    t -> boundary_condition_pointsource(applymap(param, t), center, wavenumber, amplitude)
