const PROJECT_ROOT = @__DIR__

include(joinpath(PROJECT_ROOT, "data_loader.jl"))
include(joinpath(PROJECT_ROOT, "greedy.jl"))
include(joinpath(PROJECT_ROOT, "bounds.jl"))

using .DataLoader
import .GreedyRemoveWorst: greedy_remove_worst, objective_value
using .Bounds

const SITE_ROOT = joinpath(PROJECT_ROOT, "VNSforLargeOffshoreWindFarmLayoutOptimization_SyntheticInstances", "site")
const MATRIX_ROOT = joinpath(PROJECT_ROOT, "30883835")
const DIST_ROOT = joinpath(PROJECT_ROOT, "30883835", "dist_to_substation")

const U2_VERBOSE = false

struct CableType
    cost::Float64
    capacity::Int
end

const CABLES = [
    CableType(272.0, 1),
    CableType(348.0, 2),
    CableType(480.0, 4),
    CableType(672.0, 6),
]

function cable_upper_bound(distances::Vector{Float64}, chosen::Vector{Int}, cables::Vector{CableType})
    a1 = minimum(c.cost for c in cables)
    return sum(distances[i] for i in chosen) * a1
end

function lower_bound_l1(distances::Vector{Float64}, chosen::Vector{Int}, cables::Vector{CableType})
    best = cables[argmin(c.cost / c.capacity for c in cables)]
    return sum(distances[i] for i in chosen) * (best.cost / best.capacity)
end

function get_available_zones(info)
    zones = Int[]

    for (k, v) in pairs(info.Zones)
        if String(v.Type) == "available"
            push!(zones, parse(Int, String(k)))
        end
    end

    return sort(zones)
end

function make_shifted_bound_matrix(Π::Matrix{Float64})
    n = size(Π, 1)

    offdiag_min = minimum(Π[i, j] for i in 1:n, j in 1:n if i != j)
    shift = max(0.0, -offdiag_min)

    Π_bound = copy(Π)

    for i in 1:n
        for j in 1:n
            if i != j
                Π_bound[i, j] = Π[i, j] + shift
            end
        end
    end

    return Π_bound, shift
end

function save_solution(
    zone_data::ZoneData,
    chosen_local::Vector{Int},
    z::Float64,
    U1::Float64,
    U2::Float64,
    U2_best_iter::Int,
    U2_iterations::Int,
    cable_U::Float64,
    L1::Float64
)
    chosen_global = zone_data.global_indices[chosen_local]

    out_file = "greedy_solution_u2_instance_$(zone_data.instance)_zone_$(zone_data.zone).txt"

    open(out_file, "w") do io
        println(io, "Instance: ", zone_data.instance)
        println(io, "Zone: ", zone_data.zone)
        println(io, "Target q: ", zone_data.q)
        println(io, "Selected turbines: ", length(chosen_local))
        println(io, "Objective value: ", z)
        println(io, "PDSP upper bound U1: ", U1)
        println(io, "PDSP upper bound U2: ", U2)
        println(io, "U2 iterations used: ", U2_iterations)
        println(io, "U2 best iteration: ", U2_best_iter)
        println(io, "Cable upper bound: ", cable_U)
        println(io, "Initial lower bound L1: ", L1)
        println(io, "")
        println(io, "Gap U1 - z: ", U1 - z)
        println(io, "Gap U2 - z: ", U2 - z)
        println(io, "Improvement U1 - U2: ", U1 - U2)
        println(io, "")
        println(io, "Local indices within zone:")

        for idx in chosen_local
            println(io, idx)
        end

        println(io, "")
        println(io, "Global indices in full instance:")

        for idx in chosen_global
            println(io, idx)
        end
    end

    println("Saved solution to $out_file")
end

function main()
    instances = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    # instances = ["A"]

    for instance in instances
        println("\n==============================")
        println("Loading instance $instance ...")
        println("==============================")

        data = load_instance_data(
            instance;
            site_root = SITE_ROOT,
            matrix_root = MATRIX_ROOT,
            dist_root = DIST_ROOT
        )

        zones = get_available_zones(data.info)
        println("Available zones: ", zones)

        for zone in zones
            println("\n--- Running zone $zone ---")

            zone_data = extract_zone_data(data, zone)

            println("Candidates in zone: ", length(zone_data.global_indices))
            println("Target q: ", zone_data.q)

            Π = copy(zone_data.power_submatrix)

            println("size Π = ", size(Π))
            println("min Π = ", minimum(Π), ", max Π = ", maximum(Π))
            println("symmetric? ", Π ≈ Π')

            chosen_local = greedy_remove_worst(Π, zone_data.q)
            z = objective_value(Π, chosen_local)

            Π_bound, shift = make_shifted_bound_matrix(Π)
            correction = shift * zone_data.q * (zone_data.q - 1)

            println("Off-diagonal shift used for bounds: ", shift)
            println("Correction subtracted from shifted bounds: ", correction)

            U1_bound = pdsp_u1(Π_bound, zone_data.q)

            full_iterations = size(Π_bound, 1)
            iterations = full_iterations

            U2_bound, U2_history_bound, U2_best_iter = pdsp_u2_with_history(
                Π_bound,
                zone_data.q;
                iterations = iterations,
                verbose = U2_VERBOSE
            )

            U1 = U1_bound - correction
            U2 = U2_bound - correction
            U2_history = U2_history_bound .- correction

            cable_U = cable_upper_bound(zone_data.distances, chosen_local, CABLES)
            L1 = lower_bound_l1(zone_data.distances, chosen_local, CABLES)

            println("Selected turbines: ", length(chosen_local))
            println("Objective value: ", z)
            println("PDSP upper bound U1: ", U1)
            println("PDSP upper bound U2: ", U2)
            println("U2 iterations used: ", iterations, " of ", full_iterations)
            println("U2 best iteration: ", U2_best_iter)
            println("Cable upper bound: ", cable_U)
            println("Initial lower bound L1: ", L1)
            println("Gap U1 - z: ", U1 - z)
            println("Gap U2 - z: ", U2 - z)
            println("Improvement U1 - U2: ", U1 - U2)
            println("First U2 value: ", U2_history[1])
            println("Best U2 value: ", minimum(U2_history))
            println("Last U2 value: ", U2_history[end])
            println("Unique U2 values rounded: ", length(unique(round.(U2_history, digits=8))))
            println("First 10 U2 values: ", U2_history[1:min(10, length(U2_history))])

            if U2 > U1 + 1e-8
                println("WARNING: U2 is larger than U1. This should not happen.")
            end

            save_solution(
                zone_data,
                chosen_local,
                z,
                U1,
                U2,
                U2_best_iter,
                iterations,
                cable_U,
                L1
            )
        end
    end
end

main()