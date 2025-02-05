# This script produces Table 2 of the paper 

using DiffMinPoly
using Oscar
using StructuralIdentifiability

cases = [
    rand_ode([2, 1]),
    rand_ode([2, 1, 1]),
    rand_ode([3, 1, 1]),
    rand_ode([3, 2, 2]),
    rand_ode([2, 1, 1, 1]),
    rand_ode([3, 1, 1, 1]),
    rand_ode([3, 2, 2, 2])
]
               
for c in cases
    ode = c
    x = first(sort(ode.x_vars, rev = true))
    minpoly_ord = DiffMinPoly.minpoly_order(ode, x)                

    bound_size = size(DiffMinPoly.f_min_support(ode, x, DiffMinPoly.minpoly_order(ode, x)))[1]           
    @info "Bound accordong to 'Projecting dynamical systems via a support bound' for the system "
    @info "has the size $bound_size"
end       
