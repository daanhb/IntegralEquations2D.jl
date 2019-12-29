# A collection of polygon domains

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
struct OpenLineSegmentMap{S,T} <: Parameterization{S,T}
    A   ::  SVector{2,S}       # Two endpoints
    B   ::  SVector{2,S}
end

hasparameterization(d::OpenLineSegment) = true
parameterization(d::OpenLineSegment{T}) where T = OpenLineSegmentMap{T,SVector{2,T}}(d.A, d.B)

domain(m::OpenLineSegmentMap{S}) where S = Interval{:closed,:closed,S}(0, 1)
image(m::OpenLineSegmentMap{S}) where S = OpenLineSegment{S}(m.A, m.B)

applymap(m::OpenLineSegmentMap, t) = t*m.A + (1 - t)*m.B
function gradient(m::OpenLineSegmentMap, t)
    # As gradient we take the unit vector perpendicular on AB
    # facing to the left of AB
    t = m.B - m.A
    tu = norm(t) \ t
    SVector(tu[2], tu[1])
end
