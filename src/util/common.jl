
function Base.intersect(d1::PeriodicInterval, d2::PeriodicInterval)
    @assert d1.periodicdomain == d2.periodicdomain
    n1 = numelements(d1)
    n2 = numelements(d2)
    if n1 == n2 == 1
        intersect(element(d1, 1), element(d2, 1))
    elseif n1 == 1
        intersect(element(d1, 1), element(d2, 1)) ∪ intersect(element(d1, 1), element(d2, 2))
    elseif n2 == 1
        intersect(element(d1, 1), element(d2, 1)) ∪ intersect(element(d1, 2), element(d2, 1))
    else
        L = width(d1.periodicdomain)
        domain1 = intersect(element(d1, 1), element(d2, 1)) ∪ intersect(element(d1, 1), element(d2, 2))
        domain2 = intersect(element(d1, 2), element(d2, 1)) ∪ intersect(element(d1, 2), element(d2, 2))
        PeriodicInterval(leftendpoint(domain2)..rightendpoint(domain1)+L, d1.periodicdomain)
    end
end
