include("./examples.jl")

R, (x0, x1, x2, x3, s, t, u) = polynomial_ring(QQ, ["x1", "x2", "x3", "s", "t", "u"])
F = [s^2*t + s^2*u + 4*s*t*u + 3*s*u^2 + 2*t^3 + 4*t^2*u + 2*t*u^2 + 2*u^3,
     -s^3 - 2*s^2*u - 2*s*t^2 - s*t*u + s*u^2 - 2*t*u^2 + 2*u^3,
     −s^3 − 2*s^2*t - 3*s^2*u - 3*s*t^2 - 3*s*t*u - 2*s*u^2 + 2*t^2*u - 2*t*u^2,
     s^3 + s^2*t + s^2*u - s*u^2 + t^3 + t^2*u - t*u^2 - u^3];
