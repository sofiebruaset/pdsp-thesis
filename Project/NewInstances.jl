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
