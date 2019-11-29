# @recipe function f(D::Domain)
#     title --> "Domain"
#     legend --> false
#     if hasparameterization(d)
#         param = parameterization(D)
#         supp = support(param)
#
#     else
#         error("Don't know how to plot domain ", D)
#     end
#     grid, vals
# end
