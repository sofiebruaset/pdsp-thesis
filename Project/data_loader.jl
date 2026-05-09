module DataLoader

using CodecZstd
using JSON3

export CandidatePosition,
       InstanceData,
       ZoneData,
       load_instance_data,
       extract_zone_data

struct CandidatePosition
    x::Float64
    y::Float64
    depth::Float64
    construction_cost::Float64
    zone::Int
end

struct InstanceData
    instance::String
    positions::Vector{CandidatePosition}
    distances::Vector{Float64}
    power_matrix::Matrix{Float64}
    info
end

struct ZoneData
    instance::String
    zone::Int
    global_indices::Vector{Int}
    positions::Vector{CandidatePosition}
    distances::Vector{Float64}
    power_submatrix::Matrix{Float64}
    q::Int
end

function read_space_matrix(io::IO)::Matrix{Float64}
    rows = Vector{Vector{Float64}}()
    for line in eachline(io)
        s = strip(line)
        isempty(s) && continue
        push!(rows, parse.(Float64, split(s)))
    end

    isempty(rows) && error("Matrix file is empty.")

    nrows = length(rows)
    ncols = length(rows[1])
    M = Matrix{Float64}(undef, nrows, ncols)

    for i in 1:nrows
        length(rows[i]) == ncols || error("Row $i has wrong number of columns.")
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

function load_positions(path::AbstractString)::Vector{CandidatePosition}
    positions = CandidatePosition[]
    open(path, "r") do io
        for (line_no, line) in enumerate(eachline(io))
            s = strip(line)
            isempty(s) && continue
            vals = split(s)
            length(vals) == 5 || error("availablePositions: line $line_no should have 5 columns.")
            push!(positions, CandidatePosition(
                parse(Float64, vals[1]),
                parse(Float64, vals[2]),
                parse(Float64, vals[3]),
                parse(Float64, vals[4]),
                parse(Int, vals[5])
            ))
        end
    end
    return positions
end

function load_fixed(path::AbstractString)::Vector{CandidatePosition}
    positions = CandidatePosition[]
    open(path, "r") do io
        for (line_no, line) in enumerate(eachline(io))
            s = strip(line)
            isempty(s) && continue
            vals = split(s)
            length(vals) == 5 || error("fixed_wf: line $line_no should have 5 columns.")
            push!(positions, CandidatePosition(
                parse(Float64, vals[1]),
                parse(Float64, vals[2]),
                parse(Float64, vals[3]),
                parse(Float64, vals[4]),
                parse(Int, vals[5])
            ))
        end
    end
    return positions
end

function load_distances(path::AbstractString)::Vector{Float64}
    d = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            push!(d, parse(Float64, s))
        end
    end
    return d
end

function load_info(path::AbstractString)
    return JSON3.read(read(path, String))
end

function load_instance_data(instance::String;
    site_root::String,
    matrix_root::String,
    dist_root::String)

    pos_path = joinpath(site_root, instance, "availablePositions.txt")
    fixed_path = joinpath(site_root, instance, "fixed_wf.txt")
    info_path = joinpath(site_root, instance, "info.json")
    dist_path = joinpath(dist_root, "$(instance).txt")
    matrix_path = joinpath(matrix_root, "power_matrix_instance_$(instance).txt.zst")

    candidates = load_positions(pos_path)
    fixed_pos = load_fixed(fixed_path)
    positions = vcat(candidates, fixed_pos)
    distances = vcat(load_distances(dist_path), zeros(Float64, length(fixed_pos)))
    info = load_info(info_path)
    power_matrix = load_power_matrix(matrix_path)

    length(positions) == length(distances) || error("positions and distances have different lengths.")
    size(power_matrix, 1) == length(positions) || error("power matrix and positions have inconsistent size.")
    size(power_matrix, 1) == size(power_matrix, 2) || error("power matrix must be square.")

    return InstanceData(instance, positions, distances, power_matrix, info)
end

function get_q_for_zone(info, zone::Int; target_density::Int = 4)
    zone_info = info.Zones[string(zone)]

    hasproperty(zone_info, :Type) || error("Missing Type field in info.json for zone $zone")
    String(zone_info.Type) == "available" || error("Zone $zone is not an available zone")

    hasproperty(zone_info, Symbol("Power densities MW/km2")) ||
        error("Missing 'Power densities MW/km2' in info.json for zone $zone")

    hasproperty(zone_info, Symbol("N 15MW turbines")) ||
        error("Missing 'N 15MW turbines' in info.json for zone $zone")

    densities = zone_info[Symbol("Power densities MW/km2")]
    turbines = zone_info[Symbol("N 15MW turbines")]

    idx = findfirst(x -> Int(x) == target_density, densities)
    idx === nothing && error("Power density $target_density not found for zone $zone")

    return Int(turbines[idx])
end

function extract_zone_data(data::InstanceData, zone::Int)::ZoneData
    zone_info = data.info.Zones[string(zone)]
    String(zone_info.Type) == "available" || error("Zone $zone is not an available zone")

    global_indices = [i for i in eachindex(data.positions) if data.positions[i].zone == zone]
    isempty(global_indices) && error("No positions found for zone $zone")

    positions = data.positions[global_indices]
    distances = data.distances[global_indices]

    # Start from candidate-candidate interactions inside the available zone
    Pz = copy(data.power_matrix[global_indices, global_indices])

    # Add wake contribution from fixed turbines as a diagonal shift
    fixed_indices = [i for i in eachindex(data.positions)
                     if data.positions[i].zone != zone &&
                        hasproperty(data.info.Zones[string(data.positions[i].zone)], :Type) &&
                        String(data.info.Zones[string(data.positions[i].zone)].Type) == "fixed"]

    if !isempty(fixed_indices)
        for (local_j, global_j) in enumerate(global_indices)
            fixed_effect = 0.0
            for f in fixed_indices
                fixed_effect += data.power_matrix[f, global_j] + data.power_matrix[global_j, f]
            end
            Pz[local_j, local_j] += fixed_effect
        end
    end

    q = get_q_for_zone(data.info, zone)

    return ZoneData(data.instance, zone, global_indices, positions, distances, Pz, q)
end

end