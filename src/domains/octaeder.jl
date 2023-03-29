# a regular octaeder

octa_points(::Type{T} = Float64) where T =
   SVector{3,T}[SA[1,0,0],
    SA[-1,0,0],
    SA[0,1,0],
    SA[0,-1,0],
    SA[0,0,1],
    SA[0,0,-1]
  ]

function octa_tri(::Type{T} = Float64) where T
   triangles =
      [SA[0,4,3],
       SA[3,4,1],
       SA[1,4,2],
       SA[2,4,0],
       SA[3,5,0],
       SA[1,5,3],
       SA[2,5,1],
       SA[0,5,2]]
   map(t -> t .+ 1, triangles)
end
