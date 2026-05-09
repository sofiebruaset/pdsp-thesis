include("bounds.jl")
using .Bounds

D = [
    0.0 3 7 4 10 5 7;
    3 0 9 5 5 10 6;
    7 9 0 1 3 2 4;
    4 5 1 0 1 9 1;
    10 5 3 1 0 3 2;
    5 10 2 9 3 0 3;
    7 6 4 1 2 3 0
]

q = 3

U1 = pdsp_u1(D, q)
U2, history, best_iter = pdsp_u2_with_history(D, q; iterations = size(D, 1))

println("U1 = ", U1)
println("U2 = ", U2)
println("best_iter = ", best_iter)
println("history = ", history)