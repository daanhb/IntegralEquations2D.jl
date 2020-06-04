
abstract type Parameterization{T} <: Map{T} end

gradientnorm(m::AbstractMap, t) = norm(gradient(m, t))

gradientnorm(m::DomainSets.UnitCircleMap{T}, t) where {T} = one(T)

export normal
"Evaluate the normal of the parameterized domain at the given point."
function normal(p::Parameterization, x)
    a = gradient(p, x)
    SVector(-a[2], a[1])./norm(a)
end
