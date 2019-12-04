
abstract type Parameterization{S,T} <: AbstractMap{S,T} end

jacobian(m::DomainSets.UnitCircleMap{S,T}) where {S,T} = one(S)
