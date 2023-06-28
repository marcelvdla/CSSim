# Functions for making certain plots 
using DrWatson
@quickactivate :CSSim

using PyPlot
using ChaosTools
using Random
using PyCall
using Suppressor
using LaTeXStrings

using CSSim

"""
    phase_portrait(
        T::Int, 
        params::Dict{Symbol,Any}, 
        n_points::Int = 25,
        max_density::Float64 = 4.0; 
        nrows::Int = 0,
        ncols::Int = 0,
        figsize::Tuple{Int, Int} = (8, 5))::Tuple{Figure, Matrix{PyCall.PyObject}}

Plots `params[:n]` phase portraits of `yᵢ` vs. `xᵢ` using a number `n_points`
of initial conditions and evolving the system over `T` total timesteps. Provide
`nrows` and `ncols` such that `nrows*ncols == params[:n]`.
"""
function phase_portrait(
    T::Int, 
    params::Dict{Symbol,Any}, 
    n_points::Int = 25,
    max_density::Float64 = 4.0;
    Δt = 0.01, 
    nrows::Int = 0,
    ncols::Int = 0,
    figsize::Tuple{Int, Int} = (8, 5),
    sharex::Bool = false, 
    sharey::Bool = false,
    seed::Int = 42)::Tuple{Figure, Matrix{PyCall.PyObject}}

    # Initialize data for phase portrait
    n = params[:n]
    n_states = 2
    u0s = max_density*rand(MersenneTwister(seed), n_points, n, n_states)

    @assert nrows*ncols == n "number of rows and cols must match num systems"

    # Instantiate fig+axs with the appropriate kwargs
    fig, axs = plt.subplots(
        squeeze=false, nrows=nrows, ncols=ncols, figsize=figsize,
        sharex=sharex, sharey=sharey)

    @show size(axs)

    # get the indices needed for the figure
    row_ixs, col_ixs = get_row_col_ixs_as_vectors(nrows, ncols)

    @show row_ixs
    @show col_ixs

    @assert length(row_ixs) == nrows*ncols && length(col_ixs) == nrows*ncols 
        "The number of row indices should match the number of rows*cols"

    # state variable indices
    x_ix = 1
    y_ix = 2

    # Plot phase portrait and suppress parameter warning from passing
    # mixed data types `params`
    @suppress begin 
        for initial_point in 1:n_points
            # compute trajectories 
            u0 = u0s[initial_point, :, :]
            ds = n_forest_system(u0, params)
            ds_trajectory, _ = trajectory(ds, T; Δt = Δt)

            # reshape statespace to `[timepoints, n_ecosystems, n_states]`
            ds_trajectory_tensor = reshape(
                Matrix(ds_trajectory), :, n, n_states)

            # iterate through ecosystems and plot the phase trajectories
            for ecosystem_i in 1:n
                row = row_ixs[ecosystem_i]
                col = col_ixs[ecosystem_i]
                axs[row, col].plot(
                    ds_trajectory_tensor[:, ecosystem_i, x_ix],
                    ds_trajectory_tensor[:, ecosystem_i, y_ix],
                    color = "grey")

                # Label the axes on the last point
                if initial_point == n_points
                    axs[row, col].set_xlabel(
                        LaTeXString("\$x_{$(ecosystem_i)}\$"))
                    axs[row, col].set_ylabel(
                        LaTeXString("\$y_{$(ecosystem_i)}\$"))
                end 
            end 
        end
    end 

    return fig, axs
end 

"""
    get_row_col_ixs_from_one_dim

Generate the row and column indices as flat vectors for plotting an `n`
dimensional space, handling cases for nrows==1 and ncols==1 separately.

# Examples
```julia-repl
julia> n_ecosystems = 6
julia> nrows = 3
julia> ncols = 2
julia> row_ixs, col_ixs = get_row_col_ixs_as_vectors(nrows, ncols)
([1, 1, 2, 2, 3, 3], [1, 2, 1, 2, 1, 2])
julia> for i in 1:n_ecosystems
           println(row_ixs[i], " ", col_ixs[i])
       end
1 1
1 2
2 1
2 2
3 1
3 2
```
"""
function get_row_col_ixs_as_vectors(
    nrows::Int, ncols::Int)::Tuple{Vector{Int}, Vector{Int}}
    if nrows == 1
        row_ixs = [1 for i in 1:ncols]
    else 
        row_ixs = vcat([repeat([i], ncols) for i in 1:nrows]...)
    end 

    if ncols == 1
        col_ixs = [1 for i in 1:nrows]
    else 
        col_ixs = repeat(1:ncols, nrows)
    end 

    return row_ixs, col_ixs
end 