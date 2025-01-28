# file for examples and associated convenience functions, not loaded
# by the package

using MixedFiberPolytope
using Oscar

function rand_field_elem(K)
    ch = characteristic(K)
    if !iszero(ch)
        return K(rand(1:(ch-1)))
    else
        return K(rand(1:1000))
    end
end

function square_sys_to_implicit(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    K = base_ring(R)

    S, y = polynomial_ring(QQ, "y" => 1:(ngens(R)-1))
    phi = hom(R, S, [y..., rand_field_elem(K)])

    return (phi).(F)
end

function square_sys_to_implicit(I::MPolyIdeal)
    return square_sys_to_implicit(gens(I))
end
