# This script produces Table 2 of the paper 

using MixedFiberPolytope
using Oscar
using IterTools

#-----------------------------------------------------#   

# experiments for the 3-dimentional case

R, (y, x11, x12, x13, x1, x2, x3) = polynomial_ring(QQ, ["y", "x11", "x12", "x13", "x1", "x2", "x3"])

cases = []

push!(
    cases,
    Dict(
        :name => "BlueSky",
        :polys => [(2 + 456//1000 - 10 * ((x1)^2 + (x2)^2) ) * x1 + (x2)^2 + 2 * x2 + (x3)^2,
                   - (x3)^3 - (1 + x2) * ((x2)^2 + 2 * x2 + (x3)^2) - 4 * x1 + 456//1000 * x2,
                   (1 + x2) * (x3)^2 + (x1)^2 - 357//10000],
    )
)


push!(
    cases,
    Dict(
        :name => "LV_3d",
        :polys => [x1 * (1 - 1//10 * x1 - 3//10 * x2 - 5//10 * x3) + 4//10 * (x2)^2 + 7//10 * (x3)^2,
                   x2 * (1 - 2//10 * x1 - 8//10 * x2 - 5//10 * x3 + 6//10 * (x3)^3),
                   x3 * (1 - 9//10 * x1 - 7//10 * x2 - 2//10 * x3 + 2//10 * (x3)^3)],
    )
)



for c in cases
    
    polysys = c[:polys]

    #computing the derivatives
    polys = [y - x1, x11 - polysys[1], 
             x12 - (derivative(polysys[1], x1) * x11 + derivative(polysys[1], x2) * polysys[2] + derivative(polysys[1], x3) * polysys[3])]
    push!(polys, x13 - (derivative(polys[3], x1) * x11 + derivative(polys[3], x2) * polysys[2] + derivative(polys[3], x3) * polysys[3] + derivative(polys[3], x11) * x12))
    
    #computing mixed fiber polytope
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^30)
    mf_size = length(lattice_points(mf))  
    @info "Newton polytope for $(c[:name]) has the size $mf_size and $(length(vertices(mf))) vertices"

end

#-----------------------------------------------------#

#experiments for the 4-dimentional case

R, (y, x11, x12, x13, x14, x1, x2, x3, x4) = polynomial_ring(QQ, ["y", "x11", "x12", "x13", "x14", "x1", "x2", "x3", "x4"])

cases_4dim = []

push!(
    cases_4dim,
    Dict(
        :name => "LV_4d",
        :polys => [1//10 * x1 * (6//10 + 7//10 * x1  + 3//10 * x2 - 4//10 * x3 - 3//10 * x4),
                   x2 * (2//10 + 6//10 * x1 - 3//10 * x2 - 2//10 * x3 - 1//10 * x4),
                   x3 * (1 + 5//10 * x1 + 6//10 * x2 + 8//10 * x3 + 9//10 * x4),
                   x4 * (9//10 + 7//10 * x1 + 6//10 * x2 + 2//10 * x3 + 3//10 * x4)],
    )
)

for c in cases_4dim
    
    polysys = c[:polys]

    #computing the derivatives
    polys = [y - x1, x11 - polysys[1], 
             x12 - (derivative(polysys[1], x1) * x11 + derivative(polysys[1], x2) * polysys[2] + derivative(polysys[1], x3) * polysys[3] + derivative(polysys[1], x4) * polysys[4])]
    push!(polys, x13 - (derivative(polys[3], x1) * x11 + derivative(polys[3], x2) * polysys[2] + derivative(polys[3], x3) * polysys[3] + derivative(polys[3], x4) * polysys[4]+ derivative(polys[3], x11) * x12))
    push!(
        polys, 
        x14 - ( derivative(polys[4], x1) * x11 + derivative(polys[4], x2) * polysys[2] + derivative(polys[4], x3) * polysys[3] + derivative(polys[4], x4) * polysys[4] 
                + derivative(polys[4], x11) * x12 + derivative(polys[4], x12) * x13)
    )
    
    #computing mixed fiber polytope
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^30)
    mf_size = length(lattice_points(mf))  
    @info "Newton polytope for $(c[:name]) has the size $mf_size and $(length(vertices(mf))) verts"

end

