include("./examples.jl")

R, (x0, x1, x2, x3, y1, y2, y3) = polynomial_ring(QQ, ["x0", "x1", "x2", "x3", "y1", "y2", "y3"])
F = [y1*y2 + y3 + 1,
     y1*y3 + y2 + 1,
     y1*y3 + y1 + 1,
     y1^3 + y2^5 + y3^7];
