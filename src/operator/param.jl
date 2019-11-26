
import GridArrays: GridArray

struct ParamGrid{T,N,P,A <: GridArray{T,N}} <: GridArray{T,N}
    supergrid   ::  A
    param       ::  P
end
