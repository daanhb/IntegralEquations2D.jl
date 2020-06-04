# A collection of polygon domains

export Rectangle,
    RegularPolygon,
    EquilateralTriangle,
    Square,
    Pentagon,
    Hexagon

# We first define an open line segment, this represents the side of a polygon
# given in counter clockwise order.

"An open line segment domain in `ℝ^2`."
struct OpenLineSegment{T} <: EuclideanDomain{2,T}
    A   ::  SVector{2,T}       # Two endpoints
    B   ::  SVector{2,T}

    OpenLineSegment{T}(A = SVector(-one(T), zero(T)),
            B = SVector(one(T), zero(T))) where T = new(A, B)
end

OpenLineSegment() = OpenLineSegment{Float64}()
OpenLineSegment(A::SVector{2,T}, B::SVector{2,T}) where T = OpenLineSegment{T}(A, B)
OpenLineSegment(x1::T, y1::T, x2::T, y2::T) where T = OpenLineSegment{T}(SVector(x1, y1), SVector(x2, y2))

indomain(x, ::OpenLineSegment) = error("Don't know how to evaluate indomain for an open line segment.")


"""
Parameterization of the open line segment domain, mapping from `[0,1]` to the open line segment.
"""
struct OpenLineSegmentMap{T} <: Parameterization{T}
    A   ::  SVector{2,T}       # Two endpoints
    B   ::  SVector{2,T}
end

hasparameterization(d::OpenLineSegment) = true
parameterization(d::OpenLineSegment) = OpenLineSegmentMap(d.A, d.B)

domain(m::OpenLineSegmentMap{T}) where T = zero(T)..one(T)
image(m::OpenLineSegmentMap) = OpenLineSegment(m.A, m.B)

_ols_applymap(A::SVector{2,T}, B::SVector{2,T}, t) where T =
    (1 - t)*A + t*B
applymap(m::OpenLineSegmentMap, t) = _ols_applymap(m.A, m.B, t)

function _ols_gradient(A::SVector{2,T}, B::SVector{2,T}, t) where T
    # As gradient we take the unit vector perpendicular on AB
    # facing to the left of AB
    t = B - A
    tu = norm(t) \ t
    SVector(tu[2], tu[1])
end
gradient(m::OpenLineSegmentMap, t) = _ols_gradient(A, B, t)

# Small utility function that probably is somewhere in the standard library
# Map a range to a unit range
function _map2ur(i, b, e)
    (i - b) / (e - b)
end


# A rectangle

"A rectangular domain in `ℝ^2`."
struct Rectangle{T} <: EuclideanDomain{2,T}
    ll :: SVector{2,T}   # Two endpoints
    ur :: SVector{2,T}

    function Rectangle{T}(A = SVector(-one(T), -one(T)), B = SVector(one(T), one(T))) where T
        new(SVector(min(A[1], B[1]), min(A[2], B[2])),
            SVector(max(A[1], B[1]), max(A[2], B[2])))
    end
end

Rectangle() = Rectangle{Float64}()
Rectangle(A::SVector{2,T}, B::SVector{2,T}) where T = Rectangle{T}(A, B)
Rectangle(x1::T, y1::T, x2::T, y2::T) where T = Rectangle{T}(SVector(x1, y1), SVector(x2, y2))

indomain(x::SVector{2,T}, d::Rectangle{T}) where T =
    d.ll[1] <= x[1] <= d.ur[1] && d.ll[2] <= x[2] <= d.ur[2]


"""
Parameterization of a rectangular domain, mapping from `[0,1)` to the rectangle.
"""
struct RectangleMap{T} <: Parameterization{T}
    ll :: SVector{2,T}
    ur :: SVector{2,T}
end

hasparameterization(d::Rectangle) = true
parameterization(d::Rectangle) = RectangleMap(d.ll, d.ur)

domain(m::RectangleMap{T}) where T = Interval{:closed,:open,T}(0, 1)
image(m::RectangleMap) = Rectangle(m.ll, m.ur)

function _getpart(m::RectangleMap, t)
    w, h = m.ur[1] - m.ll[1], m.ur[2] - m.ll[2]
    c = 2w + 2h

    te = w / c
    t < te && return m.ll, SVector(m.ur[1], m.ll[2]), 0, te

    ts, te = te, (w + h) / c
    t < te && return SVector(m.ur[1], m.ll[2]), m.ur, ts, te

    ts, te = te, (2w + h) / c
    t < te && return m.ur, SVector(m.ll[1], m.ur[2]), ts, te

    SVector(m.ll[1], m.ur[2]), m.ll, te, 1
end

function applymap(m::RectangleMap, t)
    s, e, st, et = _getpart(m, t)
    _ols_applymap(s, e, _map2ur(t, st, et))
end

function gradient(m::RectangleMap, t)
    s, e, st, et = _getpart(m, t)
    _ols_gradient(s, e, _map2ur(t, st, et))
end



# Regular polygons

"A regular polygon domain in `ℝ^2`."
struct RegularPolygon{T} <: EuclideanDomain{2,T}
    n      :: Int          # number of verticies
    center :: SVector{2,T} # center of polygon
    radius :: T            # the radius of the enclosing circle

    function RegularPolygon{T}(n, center=SVector(zero(T), zero(T)), radius=one(T)) where T
        if n < 3
            error("A regular polygon has 3 or more sides")
        end

        new(n, center, radius)
    end
end

RegularPolygon(n) = RegularPolygon{Float64}(n)
RegularPolygon(n::Int, center::SVector{2,T}, radius::T) where T = RegularPolygon{T}(n, center, radius)
RegularPolygon(n::Int, x::T, y::T, radius::T) where T = RegularPolygon{T}(n, SVector(x, y), radius)

function indomain(x::SVector{2,T}, d::RegularPolygon{T}) where T
    #=x -= d.center # Work relative to the origin
    Δθ = 2π / d.n
    θ0 = 3π/2 - Δθ/2=#

    # TODO
    error("Don't know how to evaluate indomain for regular polygon")
end


"""
Parameterization a regular polygon domain, mapping from `[0,1)` to the polygon.
"""
struct RegularPolygonMap{T} <: Parameterization{T}
    n      :: Int          # number of verticies
    center :: SVector{2,T} # center of the polygon
    radius :: T            # the radius of the enclosing circle
end

hasparameterization(d::RegularPolygon) = true
parameterization(d::RegularPolygon) = RegularPolygonMap(d.n, d.center, d.radius)

domain(m::RegularPolygonMap{T}) where T = Interval{:closed,:open,T}(0, 1)
image(m::RegularPolygonMap) = RegularPolygon(m.n, m.center, m.radius)

function _getpart(m::RegularPolygonMap, t)
    Δt, Δθ = 1 / m.n, 2π / m.n
    θ0 = 3π/2 - Δθ/2

    # The last segment is m.n - 1, min catches t = 1
    seg = min(floor(Int, t * m.n), m.n - 1)

    θ1 = θ0 + seg * Δθ
    y, x = sincos(θ1)
    s = m.radius * SVector(x, y) + m.center

    θ2 = θ1 + Δθ
    y, x = sincos(θ2)
    e = m.radius * SVector(x, y) + m.center

    s, e, seg * Δt, (seg + 1) * Δt
end

function applymap(m::RegularPolygonMap, t)
    s, e, st, et = _getpart(m, t)
    _ols_applymap(s, e, _map2ur(t, st, et))
end

function gradient(m::RegularPolygonMap, t)
    s, e, st, et = _getpart(m, t)
    _ols_gradient(s, e, _map2ur(t, st, et))
end



# A collection of regular polygons

EquilateralTriangle() = RegularPolygon(3)
EquilateralTriangle(center::SVector{2,T}, radius=one(T)) where T = RegularPolygon{T}(3, center, radius)
EquilateralTriangle(x::T, y::T, radius=one(T)) where T = RegularPolygon{T}(3, SVector(x, y), radius)

Square() = RegularPolygon(4)
Square(center::SVector{2,T}, radius=one(T)) where T = RegularPolygon{T}(4, center, radius)
Square(x::T, y::T, radius=one(T)) where T = RegularPolygon{T}(4, SVector(x, y), radius)

Pentagon() = RegularPolygon(5)
Pentagon(center::SVector{2,T}, radius=one(T)) where T = RegularPolygon{T}(5, center, radius)
Pentagon(x::T, y::T, radius=one(T)) where T = RegularPolygon{T}(5, SVector(x, y), radius)

Hexagon() = RegularPolygon(6)
Hexagon(center::SVector{2,T}, radius=one(T)) where T = RegularPolygon{T}(6, center, radius)
Hexagon(x::T, y::T, radius=one(T)) where T = RegularPolygon{T}(6, SVector(x, y), radius)
