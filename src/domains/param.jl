
abstract type Parameterization{S,T} <: AbstractMap{S,T} end

gradientnorm(m, t) = norm(gradient(m, t))

gradientnorm(m::DomainSets.UnitCircleMap{S,T}) where {S,T} = one(S)

struct ParametricDomain{S,T} <: DomainSets.DerivedDomain{T}
    superdomain ::  Domain{T}
    pardomain   ::  Domain{S}
    par         ::  AbstractMap{S,T}
end

ParametricDomain(d::Domain) = _ParametricDomain(d, parameterization(d))
_ParametricDomain(d, par) = ParametricDomain(d, domain(par), par)

integral(f, d::ParametricDomain; options...) = _integral(f, d, d.par; options...)

_integral(f, d::ParametricDomain, param; options...) =
    integral(t -> f(applymap(param, t))*gradientnorm(param, t), d.pardomain; options...)
