include("./examples.jl")

R, (s, t, u, v, w) = polynomial_ring(QQ, ["s", "t", "u", "v", "w"])
F = [s^3 - u^2 - t - 3*s - u + w,
     u^2 - s*w - 11,
     s^2 - 5*u - v,
     u^2 - s - v - w,
     u^2 + 7*s + t,
     v^2 + s^2 - s - t - w];
