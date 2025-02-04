module MixedFiberPolytope

using Oscar
using MixedSubdivisions
using DynamicPolynomials
using Combinatorics

export mixed_fiber_polytope, mixed_subdiv

const Support = Vector{Matrix{Int32}}
const MixCellInds = Vector{Tuple{Int, Int}}
struct MixCell
    inds::MixCellInds
    normal::Vector{Float64}
end

"""
    mixed_fiber_polytope(F::Vector{P}; implicit = false, epsinv = 50000) where {P <: Union{MPolyRingElem, Polynomial}}

If `implicit == false`, compute the mixed fiber polytope associated to
the Newton polyhedra of `F` w.r.t. the first `n-k` variables where `k
= length(F) - 1`.

If `implicit == true`, compute the mixed fiber polytope associated
to the Newton polyhedra of the system given by `t_i - F[i]` w.r.t.
new variables `t_i`. 

`epsinv` affects the likely correctness of the output, the larger it is,
the more likely the output is correct. The size of `epsinv` does not
affect the computation time but may cause overflow errors.
"""
function mixed_fiber_polytope(F::Vector{P};
                              implicit = false,
                              epsinv = 2^24) where {P <: Union{MPolyRingElem, Polynomial}}
    A = construct_support(F)
    if implicit
        @assert ambient_dim(A) == length(A) - 1 "unsuitable number of variables"
        A = implicit_support(A)
    end
    return mixed_fiber_polytope(A, epsinv = epsinv)
end

"""
    mixed_fiber_polytope(A::Vector{Matrix{Int32}})

Compute the mixed fiber polytope associated to the convex hulls of the
columns of each matrix in `A` w.r.t. the first `n-k` dimensions where
`k = length(A) - 1`.

`epsinv` affects the likely correctness of the output, the larger it is,
the more likely the output is correct. The size of `epsinv` does not
affect the computation time but may cause overflow errors.
"""
function mixed_fiber_polytope(A::Support; epsinv = 2^24)

    n = ambient_dim(A)
    k = length(A) - 1

    fiber_dict = [Dict{Int, Vector{Int}}() for _ in 1:k+1]
    pA = Matrix{Int32}[]

    # construct projected support and remember which support vectors mapped to it
    for (i, Ai) in enumerate(A)
        pAi_cols = Vector{Int32}[]
        for (j, a) in enumerate(eachcol(Ai))
            a_prj = a[n-k+1:n]
            found_im = false
            for l in 1:length(pAi_cols)
                if a_prj == pAi_cols[l]
                    found_im = true
                    add_key_or_push!(fiber_dict[i], l, j)
                    break
                end
            end
            if !found_im
                push!(pAi_cols, a_prj)
                fiber_dict[i][length(pAi_cols)] = [j]
            end
        end
        pAi = hcat(pAi_cols...)
        push!(pA, pAi)
    end

    v_orac = w -> begin
        w_rand = (qq_to_float).(w) + (rand(Float64, length(w)) ./ 10000)
        return mfp_vert(A, pA, fiber_dict, Vector(w_rand), epsinv)
    end
    
    mfp = construct_polytope(n - k, v_orac)
    return mfp
end

function mfp_vert(A::Support,
                  pA::Support,
                  fiber_dict::Vector{Dict{Int, Vector{Int}}},
                  w::Vector{Float64},
                  epsinv::Int)

    n = ambient_dim(A)
    k = length(A) - 1

    A_w = Matrix{Int32}[]

    w_ext = vcat(w, zeros(Float64, k))

    coh_weights = [Float64[] for _ in 1:k+1]
    for (i, pAi) in enumerate(pA)
        A_wi_cols = Vector{Int32}[]
        for j in eachindex(eachcol(pAi))
            preim_indices = fiber_dict[i][j]
            minval, min_idx = findmin(l -> sum(w_ext .* get_exponent(A, i, l)),
                                      preim_indices)
            preim_pt = get_exponent(A, i, preim_indices[min_idx])
            push!(A_wi_cols, preim_pt)
            push!(coh_weights[i], minval)
        end
        push!(A_w, hcat(A_wi_cols...))
    end

    # msd_weights = approximate_lifting_vector(coh_weights, epsinv)

    msd_weights = Vector{Float64}[]
    for cw in coh_weights
        push!(msd_weights, cw + rand(Float64, length(cw)) ./ epsinv)
    end

    mixed_subdiv = mixed_cells_overdet(pA, msd_weights)

    res = zeros(Int32, n-k)
    
    for (i, j, vol) in mixed_subdiv
        preim_pt = get_exponent(A_w, i, j)
        res += vol*preim_pt[1:n-k]
    end

    return res
end

# -- Helper functions for Supports -- #

function construct_support(F::Vector{P}) where {P <: MPolyRingElem}
    return (construct_support).(F)
end

function construct_support(F::Vector{P}) where {P <: Polynomial}
    return MixedSubdivisions.support(F)
end

function construct_support(f::MPolyRingElem)

    NP = convex_hull(collect(Oscar.exponents(f)))
    return hcat([Vector{Int32}((numerator).(v)) for v in vertices(NP)]...)
end

function find_minima(A::Matrix, w, tol=2^-16)

    res = [1]
    minval = sum(w .* A[:, 1])
    
    for i in 2:size(A, 2)
        evali = sum(w .* A[:, i])
        if abs(evali - minval) < tol
            push!(res, i)
        elseif evali < minval
            minval = evali
            res = [i]
        end
    end

    return res
end

function find_face(A::Vector{<:Matrix}, w)
    return [Tuple(find_minima(Ai, w)) for Ai in A]
end

function get_exponent(A::Vector{<:Matrix}, i, j)
    return A[i][:, j]
end

function ambient_dim(A::Vector{<:Matrix})
    return length(get_exponent(A, 1, 1))
end

function linear_span(A::Matrix{Int32})
    if size(A, 2) == 1
        return [zeros(Int32, size(A, 1))]
    end
    col0 = first(eachcol(A))
    return [v - col0 for v in eachcol(A)[2:end]]
end

# -- Mixed Subdivisions -- #

function vol(C::MixCell, A::Support)
    res = 1
    Cinds = C.inds
    for (i, Ci) in enumerate(Cinds)
        a, b = get_exponent(A, i, Ci[1]), get_exponent(A, i, Ci[2])
        l = Int(norm(b-a))
        res *= l
    end
    return res
end

function mixed_subdiv(A::Support, w)

    println("Computing Mixed Subdivision...")
    A_w = [vcat(Ai, transpose(wi)) for (Ai, wi) in zip(A, w)]

    cell_queue = [Tuple{Int, Int}[]]
    cell_polyhedra = [entire_sp_poly(ambient_dim(A))]

    nfeas_cache = Dict{Tuple{Int, Int, Int}, Bool}()

    res = MixCell[]
    
    k = length(A)
    while !isempty(cell_queue)

        cell, feas_pol = popfirst!(cell_queue), popfirst!(cell_polyhedra)
        l = length(cell) + 1
        if l > length(A)
            push!(res, cell)
            continue
        end

        A_wl = A_w[l]
        for (i, j) in combinations(eachindex(eachcol(A_wl)), 2)
            get(nfeas_cache, (l, i, j), false) && continue
            Pij = partial_cell_polyhedron(A_wl, i, j)
            if !is_feasible(Pij)
                nfeas_cache[(l, i, j)] = true
                continue
            end
            Q = intersect(feas_pol, Pij)
            if is_feasible(Q)
                push!(cell_queue, vcat(cell, [(i,j)]))
                push!(cell_polyhedra, Q)
            end
        end
    end

    res = MixCell[]
    for (i, C) in enumerate(res)
        aff_span = affine_span(C, A_w)
        N = nullspace(aff_span)
        size(N, 2) != 1 && continue
        v = N[:, 1]
        iszero(last(v)) && continue
        if last(v) < 0
            v = -v
        end
        v = v * 1/last(v)
        if find_face(A_w, v) == C
            push!(res, MixCell(C, v[1:end-1]))
        end
    end
        
    return res
end

function partial_cell_polyhedron(A::Matrix{Float64}, i, j)
    n = size(A, 1) - 1
    
    eq_mat = Array{Float64}(undef, 0, n)
    eq_b = Float64[]
    ineq_mat = Array{Float64}(undef, 0, n)
    ineq_b = Float64[]

    a, b = A[:, i], A[:, j]
    eq_mat = vcat(eq_mat, permutedims(a[1:n]-b[1:n]))
    push!(eq_b, b[n+1] - a[n+1])

    for (l, al) in enumerate(eachcol(A))
        l == i || l == j && continue
        ineq_mat = vcat(ineq_mat, permutedims(a[1:n] - al[1:n]))
        push!(ineq_b, al[n+1] - a[n+1])
    end

    return polyhedron((ineq_mat, ineq_b), (eq_mat, eq_b))
end

function affine_span(C::MixCellInds, A::Vector{Matrix{T}}) where T

    res = Array{T}(undef, 0, ambient_dim(A))

    for (i, Ci) in enumerate(C)
        a0, a1 = get_exponent(A, i, Ci[1]), get_exponent(A, i, Ci[2])
        res = vcat(res, permutedims(a1-a0))
    end

    return res
end

function entire_sp_poly(n)
    A = Array{Float64}(undef, 0, n)
    b = Float64[]
    return polyhedron(Float64, (A, b))
end

function mixed_cells_overdet(A::Support, w)

    k = length(A) - 1

    result = Tuple{Int, Int, Int}[]

    for i in eachindex(A)
        Ai = A[1:end .!= i]
        # is_zero_mixed_volume(Ai) && continue
        wi = w[1:end .!= i]

        # cellsi = try
        #     mixed_cells(Ai, wi)
        # catch e
        #     for Aij in Ai
        #         display(Aij)
        #     end
        #     rethrow(e)
        # end

        cellsi = mixed_subdiv(Ai, wi)

        i_supp_lifted = vcat(A[i], transpose(w[i]))

        for (j, a) in enumerate(eachcol(A[i]))
            for cell in cellsi
                nv = cell.normal
                min_indices = find_minima(i_supp_lifted, vcat(nv, [1]))
                if min_indices == [j]
                    new_cell = (i, j, vol(cell, A))
                    push!(result, new_cell)
                end
            end
        end
    end
    return result
end

# -- solve implicitization problem with given support --#

# finishing is TODO
function implicit_equation(mfp::Polyhedron, F::Vector{P}, p) where {P <: MPolyRingElem}
    supp = (Vector{Int}).(lattice_points(mfp))
    R = parent(first(F))
    n = ngens(R)
    pts = [rand(-1000:1000, n) for _ in 1:length(supp)]
    pts_im = [[qq_mod(f(pt...), p) for f in F] for pt in pts]
    mat = [prod(pt .^ mon) for pt in pts_im, mon in supp] 
    return mat
end

# -- Higgins Algorithm -- #

function construct_polytope(amb_dim::Int,
                            vert_oracle)

    start_vs = vertices(cube(amb_dim))
    vrts = unique([vert_oracle(v) for v in start_vs])
    vert_comps = length(start_vs)
    facts_confirmed = Vector{QQFieldElem}[]

    all_confirmed = false
    while !all_confirmed
        P = convex_hull(vrts)
        facts_with_nvs = [-normal_vector(fc) for fc in facets(P)]
        all_confirmed = true

        for nv in facts_with_nvs
            nv in facts_confirmed && continue
            println("vertex count: $(length(vrts))")

            new_vert = vert_oracle(nv)
            vert_comps += 1

            if !(new_vert in vrts) # check if new vertex was actually obtained
                
                push!(vrts, new_vert)
                all_confirmed = false
            else
                push!(facts_confirmed, nv)
            end
        end
    end

    println("$(vert_comps) vertex computations")

    return convex_hull(vrts)
end

# -- Other functions -- #

function qq_mod(a::QQFieldElem, p)
    F = GF(p)
    b = (Int(numerator(a) % p) * invmod(Int(denominator(a) % p), p)) % p
    return F(b)
end

function qq_to_float(a::QQFieldElem)
    return Float64(Int(numerator(a))/Int(denominator(a)))
end

function qq_to_int_vec(w::AbstractVector{QQFieldElem})
    m = lcm((denominator).(w))
    return [Int(numerator(e)) for e in m*w]
end

function add_key_or_push!(d::Dict{T, Vector{S}},
                          k::T, v::S) where {T, S}

    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = [v]
    end
end

function implicit_support(A::Support)
    k = length(A)
    n = ambient_dim(A)
    res = Matrix{Int32}[]

    for (i, Ai) in enumerate(A)
        Bi = vcat(zeros(Int32, k, size(Ai, 2)), Ai)
        ei = Int32[j == i ? 1 : 0 for j in 1:(n+k)]
        Bi = hcat(ei, Bi)
        push!(res, Bi)
    end

    return res
end

function approximate_lifting_vector(w::Vector{Vector{Float32}}, epsinv::Int)

    if all(wi -> all(iszero, wi), w)
        return (Vector{Int32}).(w)
    end

    max_entry, _ = findmax((abs).(vcat(w...)))
    multip = epsinv / max_entry

    res = Vector{Int32}[]
    for wi in w
        li = length(wi)
        rnd = rand(Int32(1):Int32(10), li) 
        push!(res, [Int32(round(wij)) for wij in multip*wi] + rnd)
    end

    return res
end

end # module MixedFiberPolytope
