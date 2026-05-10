module Bounds

export pdsp_u1, pdsp_u2, pdsp_u2_with_history

function top_k_indices(values::AbstractVector{Float64}, k::Int)::Vector{Int}
    if k < 0
        error("k must be non-negative")
    elseif k == 0
        return Int[]
    elseif k > length(values)
        error("k cannot be larger than the number of values")
    end

    return sortperm(values, rev = true)[1:k]
end

function compute_local_bounds(D::Matrix{Float64}, q::Int)
    n = size(D, 1)

    if size(D, 2) != n
        error("D must be square")
    end

    if q < 1 || q > n
        error("q must satisfy 1 <= q <= n")
    end

    dprime = zeros(Float64, n)
    y_sets = Vector{Vector{Int}}(undef, n)

    for j in 1:n
        other_indices = [i for i in 1:n if i != j]

        if q == 1
            chosen = Int[]
            value = D[j, j]
        else
            values = [D[i, j] for i in other_indices]
            pos = top_k_indices(values, q - 1)
            chosen = other_indices[pos]
            value = D[j, j] + sum(D[i, j] for i in chosen)
        end

        dprime[j] = value
        y_sets[j] = chosen
    end

    return dprime, y_sets
end

function compute_global_bound(dprime::Vector{Float64}, q::Int)
    n = length(dprime)

    if q < 1 || q > n
        error("q must satisfy 1 <= q <= n")
    end

    x_set = top_k_indices(dprime, q)
    U = sum(dprime[i] for i in x_set)

    return U, x_set
end

function pdsp_u1(D::Matrix{Float64}, q::Int)::Float64
    dprime, _ = compute_local_bounds(D, q)
    U, _ = compute_global_bound(dprime, q)
    return U
end

function pdsp_u2_with_history(
    D::Matrix{Float64},
    q::Int;
    iterations::Int = size(D, 1),
    verbose::Bool = false
)
    n = size(D, 1)

    if size(D, 2) != n
        error("D must be square")
    end

    if q < 1 || q > n
        error("q must satisfy 1 <= q <= n")
    end

    Delta = zeros(Float64, n, n)

    best_U = Inf
    best_iter = 0
    history = Float64[]

    for k in 1:iterations
        R = D .+ Delta

        dprime, y_sets = compute_local_bounds(R, q)
        U, x_set = compute_global_bound(dprime, q)

        push!(history, U)

        if U < best_U
            best_U = U
            best_iter = k
        end

        x = falses(n)
        x[x_set] .= true

        y = falses(n, n)

        for j in 1:n
            y[y_sets[j], j] .= true
        end

        newDelta = copy(Delta)

        for i in 1:n
            for j in (i+1):n
                pair_sum = D[i, j] + D[j, i]

                contribution_ij = (y[i, j] && x[j]) ? 1.0 : 0.0
                contribution_ji = (y[j, i] && x[i]) ? 1.0 : 0.0

                step_size = abs(pair_sum) / (2.0 * k)
                step = (contribution_ji - contribution_ij) * step_size

                newDelta[i, j] = Delta[i, j] + step
                newDelta[j, i] = -newDelta[i, j]
            end
        end

        Delta = newDelta

        if k <= 5
            R = D .+ Delta

            max_pair_error = maximum(abs.(
                (R[i,j] + R[j,i]) - (D[i,j] + D[j,i])
                for i in 1:n, j in 1:n if i != j
            ))

            println("iter ", k)
            println("U = ", U)
            println("max pair preservation error = ", max_pair_error)
            println("max abs Delta = ", maximum(abs.(Delta)))
        end

        if verbose
            println("U2 iteration ", k, ": ", U)
            println("Best U2 so far: ", best_U)
            println("Max abs Delta: ", maximum(abs.(Delta)))
        end
    end

    return best_U, history, best_iter
end

function pdsp_u2(
    D::Matrix{Float64},
    q::Int;
    iterations::Int = size(D, 1),
    verbose::Bool = false
)::Float64
    best_U, _, _ = pdsp_u2_with_history(
        D,
        q;
        iterations = iterations,
        verbose = verbose
    )

    return best_U
end

end