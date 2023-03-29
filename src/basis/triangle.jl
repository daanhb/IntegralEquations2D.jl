
"The map from the unit simplex to a triangle defined by three points in 3D."
function map_simplex_to_3D_triangle(p1, p2, p3)
    T = promote_type(eltype(p1),eltype(p2),eltype(p3))
    sp1 = SVector{3,T}(p1)
    sp2 = SVector{3,T}(p2)
    sp3 = SVector{3,T}(p3)
    AffineMap([sp2-sp1 sp3-sp1], sp1)
end


"Supertype of shapes functions defined on a triangle."
abstract type TriangleBasis{T} <: BasisFunctions.Dictionary{SVector{2,T},T} end

BasisFunctions.support(Φ::TriangleBasis{T}) where T = EuclideanUnitSimplex{2,T}()
BasisFunctions.measure(Φ::TriangleBasis) = lebesguemeasure(support(Φ))


"The constant shape function."
struct TriangleConstant{T} <: TriangleBasis{T}
end

TriangleConstant() = TriangleConstant{Float64}()

Base.size(Φ::TriangleConstant) = (1,)

BasisFunctions.unsafe_eval_element(Φ::TriangleConstant{T}, idx, x) where T = one(T)


"The constant and linear shape functions."
struct TriangleLinears{T} <: TriangleBasis{T}
end

TriangleLinears() = TriangleLinears{Float64}()

Base.size(Φ::TriangleLinears) = (4,)

function BasisFunctions.unsafe_eval_element(Φ::TriangleLinears{T}, idx, x) where T
    if idx == 1
        # the constant function
        one(T)
    elseif idx == 2
        # linear function equal to one in the origing
        one(T) - x[1] - x[2]
    elseif idx == 3
        # linear function equal to one in the point (1,0)
        one(T) + (x[1]-1)
    else
        # linear function equal to one in the point (0,1)
        one(T) + (x[2]-1)
    end
end
