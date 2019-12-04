# using RecipesBase

# @recipe function f(domain::IE.Kite; n=100)
#     legend --> false
#     if IE.hasparameterization(domain)
#         dpar = IE.domain(param)
#         points = EquispacedGrid(n, dpar)
#         values = applymap.(Ref(param), points)
#         values
#     else
#         error("Don't know how to plot domain ", D)
#     end
# end
