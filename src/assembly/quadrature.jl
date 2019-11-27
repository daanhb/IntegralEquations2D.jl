
using CompactTranslatesDict:
    PeriodicSubinterval

using BasisFunctions:
    unsafe_eval_element,
    unsafe_weight

# # A few routines for interoperability with pretty printing
# BasisFunctions.hasstencil(q::QuadratureStrategy) = false
# BasisFunctions.symbol(q::QuadratureStrategy) = "Quad"
# BasisFunctions.iscomposite(q::QuadratureStrategy) = false

import DomainIntegrals:
    QuadratureStrategy, QuadAdaptive,
    Q_quadgk, Q_hcubature

projectionintegral(qs, f, dict, idx, measure, sing = NoSingularity()) =
    projectionintegral(qs, f, dict, idx, measure, sing, support(dict, idx))


function projectionintegral(qs, f, dict, idx, measure, sing, domain)
    integrand = t -> f(t)*unsafe_eval_element(dict, idx, t)*unsafe_weight(measure, t)
    DomainIntegrals.integral(qs, integrand, domain, sing)
end

projectionintegral(qs, f, dict, idx, measure, sing, domain::PeriodicSubinterval) =
    sum(projectionintegral(qs, f, dict, idx, measure, sing, d) for d in elements(domain))



doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing = NoSingularity()) =
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, support(dict1, idx1), support(dict2, idx2))

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
            domain1::AbstractInterval, domain2::AbstractInterval)
    integrand = t -> f(t[1],t[2])*unsafe_eval_element(dict1, idx1, t[1])*conj(unsafe_eval_element(dict2,idx2, t[2]))*
        unsafe_weight(measure1, t[1])*unsafe_weight(measure2, t[2])
    DomainIntegrals.integral(integrand, domain1 × domain2, sing)
end

function doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing,
        domain1::PeriodicSubinterval, domain2::PeriodicSubinterval)
    sum(doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, d1, d2)
            for d1 in elements(domain1), d2 in elements(domain2))
end



"""
QuadQBF uses specialized quadrature routines that incorporate the basis function of the discretization
into their weight function. This means only evaluation of the Green's function are required.
The evaluations are equispaced and can be shared for neighbouring elements in the discretization matrix.
"""
struct QuadQBF{T} <: QuadratureStrategy
    a       ::  T       # The weight function is defined on the interval [qbf_a,qbf_b]
    b       ::  T
    x       ::  Array{T,1}  #   1d points in [qbf_a,qbf_b]
    w       ::  Array{T,1}  #   1d weights

    cbf_wr  ::  Array{T,3}  #   Weights for the regular part of the f in 2D.
    cbf_ws  ::  Array{T,3}  #   Weights for the singular part of the f in 2D.
end

QuadQBF() = QuadQBF{Float64}()

WR = zeros(5, 5, 5)
WS = zeros(5, 5, 5)

WS[:,:,1] = [
   0.000224376511000   0.001798079464820   0.000640967446610  -0.002738545656990  -0.000450592845370;
   0.004276727063400   0.047810798559800   0.047636461950690  -0.003339094770770  -0.002738545656990;
   0.008126144712220   0.105673688111370   0.132925003531680   0.047636461950690   0.000640967446610;
   0.005633676455720   0.077439254531300   0.105673688111360   0.047810798559800   0.001798079464820;
   0.000390024143850   0.005633676455720   0.008126144712220   0.004276727063400   0.000224376511000]

WS[:,:,2] = [
   0.001267080001870  -0.005967354588920  -0.013433874041170  -0.016111108889630   0.005416664932030;
  -0.006000901370660   0.025864309370310  -0.131449722949830  -0.119977127692430  -0.016111108889630;
   0.018419384832100   0.003040856819180   0.066728795114460  -0.131449722913180  -0.013433874037860;
  -0.005526798259870   0.080588550846810   0.003040856845310   0.025864309440550  -0.005967354606350;
   0.003130034264790  -0.005526798260620   0.018419384836980  -0.006000901388250   0.001267080006080]

WS[:,:,3] = [
  -0.000778533318770  -0.015009076862480   0.031568600403160  -0.033221472959820   0.012897536398830;
  -0.015009076862480  -0.102832972083900  -0.192061253278060   0.103197279120890  -0.033221472959820;
   0.031568600403160  -0.192061253278060  -0.349210902157740  -0.192061253278060   0.031568600403160;
  -0.033221472959820   0.103197279120890  -0.192061253278060  -0.102832972083900  -0.015009076862480;
   0.012897536398830  -0.033221472959820   0.031568600403160  -0.015009076862480  -0.000778533318770]

WS[:,:,4] = [
   0.001267080001870  -0.006000901370660   0.018419384832100  -0.005526798259870   0.003130034264790;
  -0.005967354588920   0.025864309370310   0.003040856819180   0.080588550846820  -0.005526798260630;
  -0.013433874041170  -0.131449722949830   0.066728795114460   0.003040856845310   0.018419384836980;
  -0.016111108889630  -0.119977127692430  -0.131449722913180   0.025864309440560  -0.006000901388250;
   0.005416664932030  -0.016111108889630  -0.013433874037860  -0.005967354606350   0.001267080006080]

WS[:,:,5] = [
   0.000224376511000   0.004276727063400   0.008126144712220   0.005633676455720   0.000390024143850;
   0.001798079464820   0.047810798559800   0.105673688111370   0.077439254531300   0.005633676455720;
   0.000640967446610   0.047636461950690   0.132925003531680   0.105673688111370   0.008126144712220;
  -0.002738545656990  -0.003339094770770   0.047636461950690   0.047810798559800   0.004276727063400;
  -0.000450592845370  -0.002738545656990   0.000640967446610   0.001798079464820   0.000224376511000]


WR[:,:,1] = [
   0.000277777777780   0.004444444444440   0.007222222222220   0.004444444444440   0.000277777777780;
   0.004444444444440   0.071111111111110   0.115555555555560   0.071111111111110   0.004444444444440;
   0.007222222222220   0.115555555555560   0.187777777777780   0.115555555555560   0.007222222222220;
   0.004444444444440   0.071111111111110   0.115555555555560   0.071111111111110   0.004444444444440;
   0.000277777777780   0.004444444444440   0.007222222222220   0.004444444444440   0.000277777777780]

WR[:,:,2] = [
  -0.000738536155200   0.007605820105820   0.002854938271610   0.011111111111110  -0.002777777777780;
   0.009938271604940   0.052460317460320   0.138633156966490   0.049761904761900   0.011111111111110;
  -0.001629188712520   0.144691358024690   0.154735449735450   0.138633156966490   0.002854938271610;
   0.011406525573190   0.047328042328040   0.144691358024690   0.052460317460320   0.007605820105820;
  -0.001715167548500   0.011406525573190  -0.001629188712520   0.009938271604940  -0.000738536155200]

WR[:,:,3] = [
  -0.001587301587300   0.012698412698410  -0.010317460317460   0.022222222222220  -0.006349206349210;
   0.012698412698410   0.044444444444440   0.161904761904760   0.025396825396820   0.022222222222220;
  -0.010317460317460   0.161904761904760   0.130158730158730   0.161904761904760  -0.010317460317460;
   0.022222222222220   0.025396825396820   0.161904761904760   0.044444444444440   0.012698412698410;
  -0.006349206349210   0.022222222222220  -0.010317460317460   0.012698412698410  -0.001587301587300]

WR[:,:,4] = [
  -0.000738536155200   0.009938271604940  -0.001629188712520   0.011406525573190  -0.001715167548500;
   0.007605820105820   0.052460317460320   0.144691358024690   0.047328042328040   0.011406525573190;
   0.002854938271600   0.138633156966490   0.154735449735450   0.144691358024690  -0.001629188712520;
   0.011111111111110   0.049761904761900   0.138633156966490   0.052460317460320   0.009938271604940;
  -0.002777777777780   0.011111111111110   0.002854938271600   0.007605820105820  -0.000738536155200]

WR[:,:,5] = [
   0.000277777777780   0.004444444444440   0.007222222222220   0.004444444444440   0.000277777777780;
   0.004444444444440   0.071111111111110   0.115555555555560   0.071111111111110   0.004444444444440;
   0.007222222222220   0.115555555555560   0.187777777777780   0.115555555555560   0.007222222222220;
   0.004444444444440   0.071111111111110   0.115555555555560   0.071111111111110   0.004444444444440;
   0.000277777777780   0.004444444444440   0.007222222222220   0.004444444444440   0.000277777777780]


QuadQBF{T}() where {T} =
    QuadQBF{T}(-T(1), T(1), [-T(1), -T(1)/2, T(0), T(1)/2, T(1)], [T(1)/60, T(16)/60, T(26)/60, T(16)/60, T(1)/60],
        WR, WS)


# Map a point x from the interval [a,b] linearly to the interval [c,d]
mapx(x, a, b, c, d) = c + (x-a)/(b-a)*(d-c)


function qbf_quadrature(f, dict, idx, measure, domain, qbf_a, qbf_b, qbf_x, qbf_w)
    a, b = extrema(domain)
    T = typeof(f(a))
    z = zero(T)
    for i = 1:length(qbf_x)
        t = mapx(qbf_x[i], qbf_a, qbf_b, a, b)
        z += qbf_w[i] * f(t) * unsafe_weight(measure, t)
    end
    z * sqrt((b-a) / (qbf_b-qbf_a))
end

function qbf_quadrature2(f, dict1, idx1, measure1, domain1, dict2, idx2, measure2, domain2, qbf_a, qbf_b, qbf_x, qbf_w)
    a1, b1 = extrema(domain1)
    a2, b2 = extrema(domain2)
    T = typeof(f(a1,a2))
    z = zero(T)
    for i = 1:length(qbf_x)
        t1 = mapx(qbf_x[i], qbf_a, qbf_b, a1, b1)
        for j = 1:length(qbf_x)
            t2 = mapx(qbf_x[j], qbf_a, qbf_b, a2, b2)
            z += qbf_w[i] * qbf_w[j] * f(t1,t2) * unsafe_weight(measure1, t1) * unsafe_weight(measure2, t2)
        end
    end
    z * sqrt((b1-a1)*(b2-a2)) / (qbf_b-qbf_a)
end


projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing::NoSingularity, domain::AbstractInterval) =
    qbf_quadrature(f, dict, idx, measure, domain, qs.a, qs.b, qs.x, qs.w)

projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing, domain::AbstractInterval) =
    projectionintegral(QuadAdaptive(), f, dict, idx, measure, sing, domain)

function projectionintegral(qs::QuadQBF, f, dict, idx, measure, sing::NoSingularity, domain::PeriodicSubinterval)
    if numelements(domain) == 1
        z = projectionintegral(qs, f, dict, idx, measure, sing, element(domain,1))
    else
        A, B = extrema(domain.periodicdomain)
        z = projectionintegral(qs, t -> f(A+mod(t-A,B-A)), dict, idx, measure, sing, domain.subdomain)
    end
    z
end


doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing::NoSingularity,
            domain1::AbstractInterval, domain2::AbstractInterval) =
        qbf_quadrature2(f, dict1, idx1, measure1, domain1, dict2, idx2, measure2, domain2, qs.a, qs.b, qs.x, qs.w)

function doubleprojection(qs::QuadQBF, f, dict1, idx1, measure1, dict2, idx2, measure2, sing::NoSingularity,
        domain1::PeriodicSubinterval, domain2::PeriodicSubinterval)

    A1, B1 = extrema(domain1.periodicdomain)
    A2, B2 = extrema(domain2.periodicdomain)
    f2 = t -> f(A1+mod(t[1]-A1,B1-A1), A2+mod(t[2]-A2,B2-A2))
    doubleprojection(qs, f, dict1, idx1, measure1, dict2, idx2, measure2, sing, domain1.subdomain, domain2.subdomain)
end



struct QuadGaussLegendre{T} <: QuadratureStrategy
    w   ::  Array{T,1}
    x   ::  Array{T,1}

    function QuadGaussLegendre{T}(n::Int) where {T}
        w, x = gausslegendre(n)
        new(x, w)
    end
end

QuadGaussLegendre(n::Int) = QuadGaussLegendre{Float64}(n::Int)
