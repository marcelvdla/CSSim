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