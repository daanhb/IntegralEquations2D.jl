
boundary_condition_wave(x, wavenumber, direction, amplitude) =
    amplitude * exp(im*wavenumber*dot(direction, x))

make_boundary_condition_planewave(wavenumber, direction, amplitude = 1) =
    (x,y) -> boundary_condition_wave(SVector(x,y), wavenumber, direction, amplitude)

make_parboundary_condition_planewave(param, wavenumber, direction, amplitude = 1) =
    t -> boundary_condition_wave(applymap(param, t), wavenumber, direction, amplitude)

boundary_condition_pointsource(x, center, wavenumber, amplitude) =
    amplitude * besselh(0, 1, wavenumber*norm(x-center))

make_boundary_condition_pointsource(center, wavenumber, amplitude = 1) =
    (x,y) -> boundary_condition_pointsource(SVector(x,y), center, wavenumber, amplitude)

make_parboundary_condition_pointsource(param, center, wavenumber, amplitude = 1) =
    t -> boundary_condition_pointsource(applymap(param, t), center, wavenumber, amplitude)
