using MixedFiberPolytope
using Oscar
using IterTools


#Randomize general polynomial of degree d over the ring R
function rand_poly(deg, vars)
    result = 0
    degs = [collect(0:deg) for v in vars]

        for m in IterTools.product(degs...)
            if sum(m) <= deg
                monom = 1 
                for i in 1:length(vars)
                    monom *= vars[i]^m[i]
                end
               result += rand(1:100) * monom
            end
        end

    return result
end

function diff_elimination_polytope(polys)
    vs = gens(parent(polys[1]))
    poly_map = [vs[1], polys[1]]
    for i in 1:(length(polys) - 1)
        newp = 0
        for (i, v) in enumerate(vs)
            newp += derivative(poly_map[end], v) * polys[i]
        end
        push!(poly_map, newp)
    end
    return poly_map
end

function my_convex_hull(points)
    pts = permutedims(hcat([Oscar.homogenize(v, 1) for v in points]...))
    ptype = Oscar._scalar_type_to_polymake(QQFieldElem)
    lin = zero_matrix(QQ, 0, size(pts, 2))
    return Oscar.Polyhedron{QQFieldElem}(Polymake.polytope.Polytope{ptype}(; VERTICES = pts, LINEALITY_SPACE=lin))
end

# R, (x1, x2, x3) = polynomial_ring(GF(65521), ["x1", "x2", "x3"])

# cases = []

# push!(
#     cases,
#     Dict(
#         :name => "Generic system [2, 1, 1]",
#         :polys => [rand_poly(2, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
#     )
# )


# for c in cases
#     polys = c[:polys]
#     print(polys)
#     P =  diff_elimination_polytope(polys)
#     print(P)
#     mf = mixed_fiber_polytope(P, implicit = true, epsinv = 50000)
#     #print(mf)
#     mf_size = length(lattice_points(mf)) 
#     print(mf_size)   
#     @info "Newton polytope for $(c[:name]) has the size $mf_size"
# end




#R, (x1111, x111, x11, y, x1, x2, x3) = polynomial_ring(GF(65521), ["x1111", "x111", "x11", "y", "x1", "x2", "x3"])

R, (y, x11, x12, x1, x2) = polynomial_ring(ZZ, ["y", "x11", "x12", "x1", "x2"])

g1 = 245*x1^2 + x1 * x2 + 2313*x2^2 + 5*x1 + 16*x2 + 1
g2 = x1 + 17*x2 + 100

f0 = y - x1

f1 = x11 - g1

f2 = x12 - ( x11 * derivative(g1, x1) + g2 * derivative(g1, x2) )

F = [f0, f1, f2]


polys_new = [sum([rand(1:500)*m for m in Oscar.monomials(p)]) for p in F] 
print(polys_new)

mf = mixed_fiber_polytope(polys_new, implicit = false, epsinv = 2^24)
mf_size = length(lattice_points(mf))  
print(mf_size)
   


# cases = []

# push!(
#     cases,
#     Dict(
#         :name => "Generic system [2, 1, 1]",
#         :polys => [rand_poly(2, [x1, x2, x3]), rand_poly(1, [x1, x2, x3]), rand_poly(1, [x1, x2, x3])],
#     )
# )

#  for c in cases

    
#     polysys = c[:polys]
#     print(polysys)
#     #computing the derivatives

#     polys = [x - x1, x11 - polysys[1], x111 - (derivative(polysys[1], x1) * x11 + derivative(polysys[1], x2) * polysys[2] + derivative(polysys[1], x3) * polysys[3])]
#     push!(polys, x1111 - (derivative(polys[3], x1) * x11 + derivative(polys[3], x2) * polysys[2] + derivative(polys[3], x3) * polysys[3] + derivative(polys[3], x11) * x111))
#     #push!(polys, x - x1)
#     print(polys)

#     #computing the generic system for polys
#     #polys_new = [sum([rand(1:500)*m for m in Oscar.monomials(p)]) for p in polys] 
#     # push!(polys_new, x1 - x)
#     # print(polys_new)  
    
    
#     mf = mixed_fiber_polytope(polys, implicit = false, epsinv = 10000)
#     mf_size = length(lattice_points(mf))  
#     print(mf_size)
#     @info "Newton polytope for $(c[:name]) has the size $mf_size"

# end







   


