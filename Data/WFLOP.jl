# run_all_greedy_removal.jl
# Run greedy removal (S=N) for all QWFLOP instances.

using DelimitedFiles
using Printf

# -----------------------------
# Settings
# -----------------------------
BASE = "instances"
M = 1e6          # penalty multiplier for constraint matrix
p = 8         # number of turbines to keep (change if your supervisor wants something else)

# -----------------------------
# Helper functions
# -----------------------------
function contribution(A::Matrix{Float64}, S::Vector{Int}, k::Int)
    c = A[k, k]
    for j in S
        if j != k
            c += A[k, j] + A[j, k]
        end
    end
    return c
end

function objective_value(A::Matrix{Float64}, S::Vector{Int})
    val = 0.0
    for i in S
        for j in S
            val += A[i, j]
        end
    end
    return val
end

function greedy_removal(Π::Matrix{Float64}, C::Matrix{Int}, p::Int; M::Float64=1e6)
    n = size(Π, 1)
    @assert size(Π, 2) == n
    @assert size(C, 1) == n && size(C, 2) == n
    @assert 1 <= p <= n

    Πtilde = Π .- M .* C
    S = collect(1:n)

    while length(S) > p
        worst_k = 0
        worst_contrib = Inf
        for k in S
            ck = contribution(Πtilde, S, k)
            if ck < worst_contrib
                worst_contrib = ck
                worst_k = k
            end
        end
        filter!(x -> x != worst_k, S)
    end

    sort!(S)
    return S, objective_value(Π, S), objective_value(Πtilde, S)
end

# -----------------------------
# Find instance folders
# -----------------------------
folders = filter(name -> startswith(name, "QWFLOP_"),
                 readdir(BASE))

sort!(folders)  # nice ordering

# -----------------------------
# Run and collect results
# -----------------------------
out_lines = String[]
push!(out_lines, "instance,n,p,objective_original,objective_penalized,selected_indices")

println("Running greedy removal for p = $p turbines")
println("Penalty multiplier M = $M")
println("------------------------------------------")

for f in folders
    path = joinpath(BASE, f)

    Πfile = joinpath(path, "cost_matrix.csv")
    Cfile = joinpath(path, "constraint_matrix.csv")

    if !(isfile(Πfile) && isfile(Cfile))
        @printf("Skipping %s (missing files)\n", f)
        continue
    end

    Π = readdlm(Πfile, ',', Float64)
    C = readdlm(Cfile, ',', Int)

    n = size(Π, 1)

    if p > n
        @printf("Skipping %s (p=%d > n=%d)\n", f, p, n)
        continue
    end

    S, objΠ, objΠtilde = greedy_removal(Π, C, p; M=M)

    @printf("%s: n=%d, p=%d, obj(Π)=%.3f\n", f, n, p, objΠ)

    # Save solution indices inside instance folder
    writedlm(joinpath(path, "greedy_removal_solution_indices.txt"), S)

    # Store CSV line
    sel = join(string.(S), " ")
    push!(out_lines, string(f, ",", n, ",", p, ",", objΠ, ",", objΠtilde, ",\"", sel, "\""))
end
