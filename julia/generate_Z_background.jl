#!/usr/bin/env julia
# filepath: /Users/strelda/Library/Mobile Documents/com~apple~CloudDocs/Files/synced Documents/phd/1 Projects/ComplexifiedZeroes/Lipkin model/generate_Z_background.jl
using LinearAlgebra, SparseArrays, DelimitedFiles, Base.Threads

# Include necessary files relative to this file.
include(joinpath(@__DIR__, "params.jl"))
include(joinpath(@__DIR__, "src", "SpinOperatorLibrary.jl"))
include(joinpath(@__DIR__, "src", "LipkinModel.jl"))
using .SpinOperatorLibrary, .LipkinModel

# Create output directory for data if it doesn't exist.
mkpath(joinpath(@__DIR__, "data"))

# Compute eigen systems.
H_initial = lipkin(lambda_initial, h_value, n_value)
H_final   = lipkin(lambda_final, h_value, n_value)

e_vals_initial, e_vecs_initial = eigen_sorted(H_initial)
e_vals_final,   e_vecs_final   = eigen_sorted(H_final)

# Construct the initial state and overlap factors.
psiI = (e_vecs_initial[:, level_initial] .+ e_vecs_initial[:, level_initial+1]) ./ sqrt(2)
# Compute all overlaps at once.
k = abs.(e_vecs_final' * psiI).^2

# Define grid based on params.jl ranges.
nbeta = background_resolution
nt = background_resolution
betas = collect(range(beta_range[1], beta_range[2], length=nbeta))
ts = collect(range(t_range[1], t_range[2], length=nt))

# Parallelized computation of the Loschmidt amplitude.
# Create a thread-local array for each thread.
local_Z = [zeros(ComplexF64, nt, nbeta) for _ in 1:nthreads()]

# Each eigenstate contributes independently.
@threads for i in eachindex(e_vals_final)
    tid = threadid()
    # Compute the outer product of two vectors for the i-th eigenstate.
    # Note: exp. returns elementwise exponentials; the multiplication below gives an nt×nbeta matrix.
    local_Z[tid] .+= k[i] .* (exp.(-im * e_vals_final[i] .* ts) *
                              exp.(-e_vals_final[i] .* betas)')
end

# Sum the contributions from each thread.
Zvals = reduce(+, local_Z)
Zvals_log = log10.(abs.(Zvals))

# Save the background data (transposing for the correct orientation)
writedlm(joinpath(@__DIR__, "data", "Z_background.txt"), Zvals_log')
println("Background data saved to data/Z_background.txt")
