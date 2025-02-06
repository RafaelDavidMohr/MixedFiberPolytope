# This script produces Table 2 of the paper 

using MixedFiberPolytope
using Oscar
using IterTools


# Randomize general polynomial of degree d over the ring R
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = 1 
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:10) * monom
            end
        end

    return result
end

#-----------------------------------------------------#   

# experiments for the 3-dimentional case

R, (y, x11, x12, x13, x1, x2, x3) = polynomial_ring(QQ, ["y", "x11", "x12", "x13", "x1", "x2", "x3"])

cases = []

push!(
    cases,
    Dict(
        :name => "BlueSky",
        :polys => [(2 + 456 - 10 * ((x1)^2 + (x2)^2) ) * x1 + (x2)^2 + 2 * x2 + (x3)^2,
        - (x3)^3 - (1 + x2) * ((x2)^2 + 2 * x2 + (x3)^2) - 4 * x1 + 456 * x2,
        (1 + x2) * (x3)^2 + (x1)^2 - 357],
    )
)



push!(
    cases,
    Dict(
        :name => "LV_3d",
        :polys => [x1 * (1 - 10 * x1 - 3 * x2 - 5 * x3) + 4 * (x2)^2 + 7 * (x3)^2,
        x2 * (1 - 2 * x1 - 8 * x2 - 5 * x3 + 6 * (x3)^3),
        x3 * (1 - 9 * x1 - 7 * x2 - 2 * x3 + 2 * (x3)^3)],
    )
)



 for c in cases
    
    polysys = c[:polys]

    #computing the derivatives
    polys = [y - x1, x11 - polysys[1], 
        x12 - (derivative(polysys[1], x1) * x11 + derivative(polysys[1], x2) * polysys[2] + derivative(polysys[1], x3) * polysys[3])]
    push!(polys, x13 - (derivative(polys[3], x1) * x11 + derivative(polys[3], x2) * polysys[2] + derivative(polys[3], x3) * polysys[3] + derivative(polys[3], x11) * x12))
   
    #computing mixed fiber polytope
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^29)
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
        :polys => [x1 * (12 - 8 * x1  + 8 * x2 - 7 * x3 - 10 * x4),
         x2 * (13 + 7 * x1 - 5 * x2 - 7 * x3 - 9 * x4),
          x3 * (-19 + x3 + 8 * x1 + 5 * x2 + x4),
           x4 * (-9 + x4 + 17 * x1 + 6 * x2 + 9 * x3)],
    )
)

# push!(
#     cases_4dim,
#     Dict(
#         :name => "Generic system [3, 1, 1, 1]",
#         :polys => [rand_poly(3, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4])],
#     )
# )

# push!(
#     cases_4dim,
#     Dict(
#         :name => "Generic system [3, 2, 2, 2]",
#         :polys => [rand_poly(3, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4])],
#     )
# )


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
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^24)
    mf_size = length(lattice_points(mf))  
    @info "Newton polytope for $(c[:name]) has the size $mf_size"

end

