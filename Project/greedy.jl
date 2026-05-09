module GreedyRemoveWorst

export greedy_remove_worst, objective_value, pdsp_u1

function initial_contributions(Π::Matrix{Float64})::Vector{Float64}
    row_sums = vec(sum(Π, dims=2))
    col_sums = vec(sum(Π, dims=1))

    n = size(Π, 1)
    contrib = zeros(Float64, n)

    for i in 1:n
        contrib[i] = row_sums[i] + col_sums[i] - Π[i, i]
    end

    return contrib
end

function greedy_remove_worst(Π::Matrix{Float64}, q::Int)::Vector{Int}
    n = size(Π, 1)
    1 <= q <= n || error("q must satisfy 1 <= q <= n")

    selected = trues(n)
    contrib = initial_contributions(Π)
    remaining = n

    while remaining > q
        worst_i = 0
        worst_val = Inf

        for i in 1:n
            if selected[i]
                val = contrib[i]
                if val < worst_val || (val == worst_val && (worst_i == 0 || i < worst_i))
                    worst_val = val
                    worst_i = i
                end
            end
        end

        selected[worst_i] = false
        remaining -= 1

        for j in 1:n
            if selected[j]
                contrib[j] -= Π[j, worst_i] + Π[worst_i, j]
            end
        end
    end

    return findall(selected)
end

function objective_value(Π::Matrix{Float64}, chosen::Vector{Int})::Float64
    return sum(Π[chosen, chosen])
end

function pdsp_u1(Π::Matrix{Float64}, q::Int)::Float64
    n = size(Π, 1)

    1 <= q <= n || error("q must satisfy 1 <= q <= n")

    dprime = zeros(Float64, n)

    for j in 1:n
        values = Float64[]

        for i in 1:n
            if i != j
                push!(values, Π[i, j])
            end
        end

        sort!(values, rev = true)

        best_sum = sum(values[1:q-1])

        dprime[j] = Π[j, j] + best_sum
    end

    sort!(dprime, rev = true)

    return sum(dprime[1:q])
end

end