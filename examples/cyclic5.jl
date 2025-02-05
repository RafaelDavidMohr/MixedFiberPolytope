include("./examples.jl")

R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x$i" for i in 1:5])

G = [x1 + x2 + x3 + x4 + x5,
     x1*x2 + x1*x5 + x2*x3 + x3*x4 + x4*x5,
     x1*x2*x3 + x1*x2*x5 + x1*x4*x5 + x2*x3*x4 + x3*x4*x5,
     x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5,
     x1*x2*x3*x4*x5 - 1]
F = square_sys_to_implicit(G);
