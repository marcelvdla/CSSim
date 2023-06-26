using DrWatson
@quickactivate :CSSim

using PyPlot
using ChaosTools
using Random
using PyCall
using Suppressor
using LaTeXStrings

include(joinpath(srcdir(), "n_forest.jl"))

"""
    phase_portrait(
        T::Int, 
        params::Dict{Symbol,Any}, 
        n_points::Int = 25,
        max_density::Float64 = 4.0; 
        plot_fixed_points::Bool = false,
        nrows::Int = 0,
        ncols::Int = 0,
        figsize::Tuple{Int, Int} = (8, 5))::Tuple{Figure, Matrix{PyCall.PyObject}}

Plots `params[:n]` phase portraits of `yᵢ` vs. `xᵢ`

Optionally `plot_fixed_points` on phase portraits.
"""
function phase_portrait(
    T::Int, 
    params::Dict{Symbol,Any}, 
    n_points::Int = 25,
    max_density::Float64 = 4.0; 
    plot_fixed_points::Bool = false,
    nrows::Int = 0,
    ncols::Int = 0,
    figsize::Tuple{Int, Int} = (8, 5))::Tuple{Figure, Matrix{PyCall.PyObject}}

    # Initialize data for phase portrait
    n = params[:n]
    n_states = 2
    u0s = max_density*rand(MersenneTwister(42), n_points, n, n_states)

    @assert nrows*ncols == n "number of rows and cols must match num systems"

    # Instantiate fig+axs with the appropriate kwargs
    fig, axs = plt.subplots(
        squeeze=false, nrows=nrows, ncols=ncols, figsize=figsize)

    @show size(axs)

    # Generate the row and column indices as flat vectors for the plotting 
    # ... handles case for nrows==1 and ncols==1 separately
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
            u0 = u0s[initial_point, :, :]
            ds = n_forest_system(u0, params)
            ds_trajectory, _ = trajectory(ds, T)

            # reshape statespace to `[timepoints, n_ecosystems, n_states]`
            ds_trajectory_tensor = reshape(
                Matrix(ds_trajectory), :, n, n_states)

            # iterate through ecosystems and plot the phase trajectories
            # TODO: move to function for performance...
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

        # TODO: Compute fixed points??
        
    end 

    return fig, axs
end 
