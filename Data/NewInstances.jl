base_dir = "30883835"

try
    import CodecZstd
catch
    println("CodecZstd is not installed. Run:")
    println("  julia -e 'import Pkg; Pkg.add(\"CodecZstd\")'")
    error("Cannot read .zst files without CodecZstd")
end

function read_space_matrix(io::IO)
    rows = Vector{Vector{Float64}}()
    for line in eachline(io)
        line = strip(line)
        isempty(line) && continue
        push!(rows, parse.(Float64, split(line)))
    end
    if isempty(rows)
        return zeros(Float64, 0, 0)
    end

    nrows = length(rows)
    ncols = length(rows[1])
    for (i, row) in enumerate(rows)
        if length(row) != ncols
            error("Line $i has length $(length(row)), expected $ncols")
        end
    end

    matrix = Matrix{Float64}(undef, nrows, ncols)
    for i in 1:nrows
        matrix[i, :] = rows[i]
    end
    return matrix
end

function load_matrix_file(path::AbstractString)
    println("Loading power matrix from: $path")
    open(path) do f
        if endswith(path, ".zst")
            stream = CodecZstd.ZstdDecompressorStream(f)
            return read_space_matrix(stream)
        else
            return read_space_matrix(f)
        end
    end
end

function load_all_power_matrices(base_dir::AbstractString)
    files = sort(filter(f -> startswith(f, "power_matrix_instance_") && endswith(f, ".txt.zst"), readdir(base_dir)))
    if isempty(files)
        error("No power matrix files found in $base_dir")
    end

    matrices = Dict{String, Matrix{Float64}}()
    for fname in files
        key = replace(fname, r"^power_matrix_instance_([A-Z])\.txt\.zst" => s"\1")
        path = joinpath(base_dir, fname)
        matrices[key] = load_matrix_file(path)
        println("Loaded instance $key size = ", size(matrices[key]))
    end
    return matrices
end

matrices = load_all_power_matrices(base_dir)
println("Loaded $(length(matrices)) power matrices:")
for key in sort(collect(keys(matrices)))
    println(" - instance $key: size = ", size(matrices[key]))
end

using CodecZstd

function read_space_matrix(io::IO)
    rows = Vector{Vector{Float64}}()
    for line in eachline(io)
        s = strip(line)
        isempty(s) && continue
        push!(rows, parse.(Float64, split(s)))
    end

    nrows = length(rows)
    ncols = length(rows[1])
    M = Matrix{Float64}(undef, nrows, ncols)

    for i in 1:nrows
        length(rows[i]) == ncols || error("Row $i has wrong number of columns")
        M[i, :] = rows[i]
    end
    return M
end

function load_power_matrix(path::AbstractString)::Matrix{Float64}
    open(path, "r") do f
        if endswith(path, ".zst")
            stream = CodecZstd.ZstdDecompressorStream(f)
            return read_space_matrix(stream)
        else
            return read_space_matrix(f)
        end
    end
end

function build_score_matrix(P::Matrix{Float64}; cable_estimate::Union{Nothing,Vector{Float64}}=nothing)
    Π = copy(P)
    if cable_estimate !== nothing
        n = size(Π, 1)
        length(cable_estimate) == n || error("cable_estimate has wrong length")
        for i in 1:n
            Π[i, i] -= cable_estimate[i]
        end
    end
    return Π
end

function initial_contributions(Π::Matrix{Float64})::Vector{Float64}
    n = size(Π, 1)
    contrib = zeros(Float64, n)
    row_sums = vec(sum(Π, dims=2))
    col_sums = vec(sum(Π, dims=1))
    for i in 1:n
        contrib[i] = row_sums[i] + col_sums[i] - Π[i, i]
    end
    return contrib
end

function greedy_remove_worst(Π::Matrix{Float64}, q::Int)
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
                if val < worst_val || (val == worst_val && i < worst_i)
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