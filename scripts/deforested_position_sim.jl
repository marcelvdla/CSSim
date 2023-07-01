# Run a position based simulation
using DataFrames: DataFrame
using DynamicalSystems: trajectory
using ProgressBars: ProgressBar
using Random: MersenneTwister
using Suppressor
using CSSim

function deforested_position_sim(
    n_sims::Int, 
    n_ecosystems::Int = 10; 
    min_y_density::AbstractFloat = 4.0,
    max_x_density::AbstractFloat = 4.0)

    # run simulation 
    print("Running deforestion simulation:")
    sim_results = simulation(
        n_sims, 
        n_ecosystems; 
        max_x_density=max_x_density, 
        min_y_density=min_y_density)

    # save simulation results 
    println("Simulation complete!")
    return DataFrame(sim_results)
end

"""
    simulation(
        n_sims::Int, 
        n_ecosystems::Int, 
        T::Int = 100, 
        T_deforest::Int = 40;
        min_y_density::Float64,
        max_x_density::Float64,
        D::Float64 = 900.0,
        y_ix::Int = 2)::Dict{Symbol, Vector}

Run deforestation simulation for `n_sims` for `n_ecosystems` by evolving
the system for `T` timesteps and deforesting `N_deforested` ecosystems in 
some random integer range in `T_deforest`. The `max_denesity` is the density 
used to initialize the old tree species densities (i.e., `yᵢ`).
"""
function simulation(
    n_sims::Int, 
    n_ecosystems::Int, 
    T::Int = 100, 
    T_deforest::Int = 40;
    min_y_density::Float64,
    max_x_density::Float64,
    D::Float64 = 900.0,
    y_ix::Int = 2)::Dict{Symbol, Vector}

    # set default parameters 
    params =  Dict{Symbol,Any}(
        :ρ  => 4.2, 
        :f  => 1.0,
        :α₀ => -1.0, 
        :w₀ =>  1.0,
        :a₁ => 1.0, 
        :h  => 2.0, 
        :a₂ => 0.0, 
        :d  => repeat([D / (n_ecosystems-1)], n_ecosystems-1), 
        :l  => 600.0, 
        :P₀ => 1.05, 
        :β₁ => 0, 
        :β₂ => 0.15,
        :n  => n_ecosystems,
        :ecosystems_to_deforest => undef,
    )

    sim_results = Dict{Symbol, Vector}(
        :num_ecosystems_killed => [],
        :deforested_position_1 => [],
        :deforested_position_2 => [],
    )

    # For n-simulations, deforest two such forests for all possible positions
    # in the network of n-ecosystems 
    for sim in ProgressBar(1:n_sims)
        for deforested_position_1 in 1:n_ecosystems
            for deforested_position_2 in deforested_position_1+1:n_ecosystems
                # Set experimental parameters
                seed = sim 
                random_tstars_place_holder = []
                params[:ecosystems_to_deforest] = ecosystems_times_to_deforest(
                        T_deforest,
                        [deforested_position_1, deforested_position_2],
                        random_tstars_place_holder,
                        n_ecosystems,
                        seed
                    )
                
                # Initial conditions
                rng = MersenneTwister(sim)
                y₀ = min_y_density .+ rand(rng, n_ecosystems)
                x₀ = max_x_density*rand(rng, n_ecosystems)
                u₀ = hcat(x₀, y₀)

                # Evolve the system through time and collect the results 
                # as a [integration_timesteps, n_ecosystems, n_states] tensor 
                time_ecosystem_state_tensor = experiment(u₀, params, T)

                # Compute the number of killed ecosystems 
                # this particular configuration of deforested positions 
                ecosystem_ylast = time_ecosystem_state_tensor[
                    end, :, y_ix]

                # A bit vector of length [n_ecosystems] where each element 
                # is 0 for a killed ecosytem or 1 if it survived 
                ecosystems_killed_vector = iskilled.(ecosystem_ylast)

                # Compute the number of ecosystems killed during this run 
                # of the simulatio
                ecosystems_killed = sum(ecosystems_killed_vector)

                # update the results dictionary
                push!(sim_results[:num_ecosystems_killed], ecosystems_killed)
                push!(
                    sim_results[:deforested_position_1], deforested_position_1)
                push!(
                    sim_results[:deforested_position_2], deforested_position_2)
            end 
        end 
    end

    return sim_results
end 

function experiment(u0, params::Dict{Symbol, Any}, T::Int)
    # evolve the n-forest system over time 
    @suppress begin
        ds = n_forest_system(u0, params)
        X, t = trajectory(ds, T)

        # define dimensions for reshaping 
        n_states = 2
        time = : ;
        n_ecosystems = params[:n]

        # reshape the data to more interpretable dimensions 
        time_ecosystem_state_tensor = reshape(
            Matrix(X), time, n_ecosystems, n_states)

        return time_ecosystem_state_tensor
    end 
end 
