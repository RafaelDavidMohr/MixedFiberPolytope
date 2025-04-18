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

R, (y, x11, x12, x13, x1, x2, x3) = polynomial_ring(ZZ, ["y", "x11", "x12", "x13", "x1", "x2", "x3"])

cases = []

push!(
    cases,
    Dict(
        :name => "Generic system [2, 1, 1]",
        :polys => [rand_poly(2, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
    )
)


push!(
    cases,
    Dict(
        :name => "Generic system [3, 1, 1]",
        :polys => [rand_poly(3, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
    )
)

push!(
    cases,
    Dict(
        :name => "Generic system [3, 2, 2]",
        :polys => [rand_poly(3, [x1, x2, x3]), rand_poly(2, [x1, x2, x3]), rand_poly(2, [x1, x2, x3])],
    )
)



 for c in cases
    
    polysys = c[:polys]

    #computing the derivatives
    polys = [y - x1, x11 - polysys[1], 
        x12 - (derivative(polysys[1], x1) * x11 + derivative(polysys[1], x2) * polysys[2] + derivative(polysys[1], x3) * polysys[3])]
    push!(polys, x13 - (derivative(polys[3], x1) * x11 + derivative(polys[3], x2) * polysys[2] + derivative(polys[3], x3) * polysys[3] + derivative(polys[3], x11) * x12))
   
    #computing mixed fiber polytope
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^24)
    mf_size = length(lattice_points(mf))  
    @info "Newton polytope for $(c[:name]) has the size $mf_size"

end

#-----------------------------------------------------#

# experiments for the 4-dimentional case

R, (y, x11, x12, x13, x14, x1, x2, x3, x4) = polynomial_ring(ZZ, ["y", "x11", "x12", "x13", "x14", "x1", "x2", "x3", "x4"])

cases_4dim = []

push!(
    cases_4dim,
    Dict(
        :name => "Generic system [2, 1, 1, 1]",
        :polys => [rand_poly(2, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4])],
    )
)

push!(
    cases_4dim,
    Dict(
        :name => "Generic system [3, 1, 1, 1]",
        :polys => [rand_poly(3, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4]), rand_poly(1, [x1, x2, x3, x4])],
    )
)

push!(
    cases_4dim,
    Dict(
        :name => "Generic system [3, 2, 2, 2]",
        :polys => [rand_poly(3, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4]), rand_poly(2, [x1, x2, x3, x4])],
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
    mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 2^24)
    mf_size = length(lattice_points(mf))  
    @info "Newton polytope for $(c[:name]) has the size $mf_size"

end

