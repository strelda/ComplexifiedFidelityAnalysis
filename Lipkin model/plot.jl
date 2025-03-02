#!/usr/bin/env julia
using LinearAlgebra, SparseArrays, DelimitedFiles, Plots

plotly()

include(joinpath(@__DIR__, "params.jl"))

# Load precomputed background values of Z.
Zvals_log = readdlm(joinpath(@__DIR__, "data", "Z_background.txt"))

# Recreate the grid using the same parameters.
nbeta = size(Zvals_log, 2)
nt = size(Zvals_log, 1)
betas = range(beta_range[1], beta_range[2], length=nbeta)
ts = range(t_range[1], t_range[2], length=nt)

# Load the approximate roots from file.
approx_roots_data = readdlm(joinpath(@__DIR__, "data", "approx_roots.txt"))
approx_roots = [complex(row[1], row[2]) for row in eachrow(approx_roots_data)]

# Plot the heatmap and overlay the approximate zero locations.
# Use a layout that forces the x-axis scale to be anchored to the y-axis.
plt = heatmap(betas, ts, Zvals_log,
              color=:viridis, xlabel="β", ylabel="t",
              xlims=(beta_range[1], beta_range[2]),
              ylims=(t_range[1], t_range[2]),
              size=(800,800))

scatter!(plt, real.(approx_roots), imag.(approx_roots),
         marker=:circle, markersize=1, color=:red, label="", markerstrokecolor=:red)

# Construct the filename with parameter information
filename = "loschmidt_plot_lambda_i=$(lambda_initial)_lambda_f=$(lambda_final)_h=$(h_value)_n=$(n_value)_recursion=$(recursion_steps)"

# Save the plot to PDF and interactive HTML.
savefig(plt, joinpath(@__DIR__, "../img_temporary", filename * ".pdf"))
# savefig(plt, joinpath(@__DIR__, "img", filename * ".html"))
println("Plots saved to the ../img_temporary folder.")