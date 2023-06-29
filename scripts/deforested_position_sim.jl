#=
# Pseudo code for running position based simulation
data = []
for sim in 1:n_sims
    for p1 in 1:10
        for p2 in 1:10
            seed = sim
            y0 = 
            x0 = rand(seed)
            u0 = [x0, y0] # not totally right 
            ecosystems_to_deforest(..., seed)
            X, t = trajectory(ds)
            push!(data, X)
        end 
    end 
end 
=#
using Random
using ProgressBars
using CSSim


"""
    simulation(
        n_sims::Int, 
        n_ecosystems::Int, 
        N_deforested::Int, 
        T::int, 
        T_deforest::Int)

Run deforestation simulation for `n_sims` for `n` forest ecosystems by evolving
the system for `T` timesteps and deforesting `N` ecosystems in some random 
integer range in `T_deforest`. The `max_denesity` is the density used 
to initialize the old tree species densities (i.e., `yᵢ`)
"""
function simulation(
    n_sims::Int, 
    n_ecosystems::Int, 
    N_deforested::Int, 
    T::int, 
    T_deforest::Int;
    outfile::String,
    max_y_density::float = 2.25,
    max_x_density::float = 4.0)

    n_states = 2
    trajectories_for_sims = []

    println("Running deforestation simulation:") 
    for sim in ProgressBar(1:n_sims)
        for deforested_position_1 in 1:n_ecosystems
            for deforested_position_2 in deforested_position_1+1:n_ecosystems
                params::Dict{Symbol, Any} = Dict(
                    # TODO
                )
                y₀ = max_y_density*rand(MersenneTwister(sim), n_ecosystems)
                x₀ = max_x_density*rand(MersenneTwister(sim), n_ecosystems)
                u₀ = hcat(x₀, y₀)
                trajectories_tensor = experiment(u₀, params, T)
                push!(trajectories_for_sims, trajectories_tensor)
            end 
        end 
    end 

    println("Writing to $(outfile)")
    # ....
end 

function experiment(u0, params::Dict{Symbol, Any}, T::Int)

end 

function main()

end 

main()