include("./examples.jl")

R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["x$i" for i in 1:6])
G = [x1+2*x2+2*x3+2*x4+2*x5+2*x6-1,
     x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2-x1,
     2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6-x2,
     x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6-x3,
     2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6-x4,
     x3^2+2*x2*x4+2*x1*x5+2*x2*x6-x5]

F = square_sys_to_implicit(G);
