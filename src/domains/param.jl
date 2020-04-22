
abstract type Parameterization{S,T} <: AbstractMap{S,T} end

jacobian(m::DomainSets.UnitCircleMap{S,T}) where {S,T} = one(S)

export normal
"Evaluate the normal of the parameterized domain at the given point."
function normal(p::Parameterization, x)
    a = gradient(p, x)
    SVector(-a[2], a[1])./norm(a)
end
