include("data_loader.jl")
include("greedy.jl")

using .DataLoader
using .GreedyRemoveWorst

const SITE_ROOT = "VNSforLargeOffshoreWindFarmLayoutOptimization_SyntheticInstances/site"
const MATRIX_ROOT = "30883835"
const DIST_ROOT = "30883835/dist_to_substation"


function get_available_zones(info)
    zones = Int[]
    for (k, v) in pairs(info.Zones)
        if String(v.Type) == "available"
            push!(zones, parse(Int, String(k)))
        end
    end
    return sort(zones)
end

function save_solution(zone_data::ZoneData, chosen_local::Vector{Int}, z::Float64)
    chosen_global = zone_data.global_indices[chosen_local]

    out_file = "greedy_solution_instance_$(zone_data.instance)_zone_$(zone_data.zone).txt"
    open(out_file, "w") do io
        println(io, "Instance: ", zone_data.instance)
        println(io, "Zone: ", zone_data.zone)
        println(io, "Target q: ", zone_data.q)
        println(io, "Selected turbines: ", length(chosen_local))
        println(io, "Objective value: ", z)
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

function parse_args()
    instance = length(ARGS) >= 1 ? ARGS[1] : "J"
    zone = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : nothing
    return instance, zone
end

function main()
    #instances = ["A","B","C","D","E","F","G","H","I","J"]
    instances = ["A"]
    for instance in instances
        println("\n==============================")
        println("Loading instance $instance ...")
        println("==============================")

        data = load_instance_data(instance;
            site_root=SITE_ROOT,
            matrix_root=MATRIX_ROOT,
            dist_root=DIST_ROOT
        )

        zones = get_available_zones(data.info)
        println("Available zones: ", zones)

        for zone in zones
            println("\n--- Running zone $zone ---")

            zone_data = extract_zone_data(data, zone)

            println("Candidates in zone: ", length(zone_data.global_indices))
            println("Target q: ", zone_data.q)

            Π = zone_data.power_submatrix

            chosen_local = greedy_remove_worst(Π, zone_data.q)
            z = objective_value(Π, chosen_local)

            println("Selected turbines: ", length(chosen_local))
            println("Objective value: ", z)

            save_solution(zone_data, chosen_local, z)
        end
    end
end


main()
