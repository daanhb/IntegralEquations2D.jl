
# A collection of smooth domains

## A kite-shaped domain

export Kite
"A kite-shaped domain in `ℝ^2`."
struct Kite{T} <: EuclideanDomain{2,T}
    a   ::  Int       # a parameter
end

Kite() = Kite{Float64}()
Kite{T}() where {T} = Kite{T}(2)
Kite(a::Integer) = Kite{float(typeof(a))}(a)

indomain(x, ::Kite) = error("Don't know how to evaluate indomain for a kite.")


using BasisFunctions: subeltype

"""
Smooth parameterization of the kite-shaped domain, mapping from `[0,1)]` to the kite.
"""
struct KiteMap{T} <: Parameterization{T}
    a   ::  Int       # a parameter
end

hasparameterization(d::Kite) = true
parameterization(d::Kite{T}) where T = KiteMap{T}(d.a)

domain(d::KiteMap{T}) where T = Interval{:closed,:open,T}(0, 1)
image(m::KiteMap{T}) where T = Kite{T}(m.a)

function applymap(m::KiteMap, t)
    b = 2*convert(domaintype(m), pi)
    SVector(-sin(b*t)-cos(m.a*b*t), cos(b*t))
end

function gradient(m::KiteMap, t)
    b = 2*convert(domaintype(m), pi)
    SVector(-b*cos(b*t)+m.a*b*sin(m.a*b*t), -b*sin(b*t))
end

function gradient_derivative(m::KiteMap, t)
    b = 2*convert(domaintype(m), pi)
    SVector(b^2*sin(b*t)+(m.a*b)^2*cos(m.a*b*t), -b^2*cos(b*t))
end



## An ellipse

export Ellipse
"An ellipse-shaped domain"
struct Ellipse{T} <: EuclideanDomain{2,T}
    center  ::  SVector{2,T}        # center of the ellipse
    radius  ::  SVector{2,T}        # (radius1,radius2)

    Ellipse{T}(center = SVector(zero(T),zero(T)), radius = SVector(one(T),one(T))) where {T} = new(center, radius)
end

Ellipse() = Ellipse{Float64}()

Ellipse(center::SVector{2,T}, radius::SVector{2,T}) where {T} = Ellipse{T}(center, radius)

Ellipse(x, y, radius1, radius2) = Ellipse(SVector(x, y), SVector(radius1, radius2))

indomain(x, d::Ellipse) = ((x[1]-d.center[1])/d.radius[1])^2 + ((x[2]-d.center[2])/d.radius[2])^2 ≈ 1


"""
Smooth parameterization of an ellipse, mapping from `[0,1)]`.
"""
struct EllipseMap{T} <: Parameterization{T}
    center  ::  SVector{2,T}        # center of the ellipse
    radius  ::  SVector{2,T}        # (radius1,radius2)
end

hasparameterization(d::Ellipse) = true
parameterization(d::Ellipse) = EllipseMap(d.center, d.radius)

domain(m::EllipseMap{T}) where T = Interval{:closed,:open,T}(0, 1)
image(m::EllipseMap) = Ellipse(m.center, m.radius)

function applymap(m::EllipseMap, t)
    b = 2*convert(domaintype(m), pi)
    m.center + SVector(m.radius[1]*cos(b*t), m.radius[2]*sin(b*t))
end

function gradient(m::EllipseMap, t)
    b = 2*convert(domaintype(m), pi)
    SVector(-m.radius[1]*b*sin(b*t),  m.radius[2]*b*cos(b*t))
end

function gradient_derivative(m::EllipseMap, t)
    b = 2*convert(domaintype(m), pi)
    SVector(m.radius[1]*b^2*cos(b*t),  m.radius[2]*b^2*sin(b*t))
end


# # function normal{T}(b::Ellipse{T}, t)
# #     n = Vec2{T}(b.radius2 * cos(2*T(pi)*t), b.radius1 * sin(2*T(pi)*t))
# #     n / norm(n)
# # end
#
# gradient_d2{T}(b::Ellipse{T}, t) = Vec2{T}(-b.radius1 * (2 * T(pi))^2 * cos(2*T(pi)*t),  -b.radius2 * (2 * T(pi))^2 * sin(2*T(pi)*t))
# gradient_d3{T}(b::Ellipse{T}, t) = Vec2{T}( b.radius1 * (2 * T(pi))^3 * sin(2*T(pi)*t),  -b.radius2 * (2 * T(pi))^3 * cos(2*T(pi)*t))
