#!/usr/bin/env julia
# filepath: /Users/strelda/Library/Mobile Documents/com~apple~CloudDocs/Files/synced Documents/phd/1 Projects/ComplexifiedZeroes/Lipkin model/find_zeroes_cont.jl

using LinearAlgebra, SparseArrays, DelimitedFiles, Base.Threads, QuadGK, FiniteDiff

# -------------------------------
# Timing & Logging
# -------------------------------
start_time = time()

global detailed_log_file = joinpath(@__DIR__, "data", "detailed_log.txt")
mkpath(joinpath(@__DIR__, "data"))
open(detailed_log_file, "w") do f
    println(f, "=== Detailed Log Start ===")
end

"""
    log_info(msg::String)

Appends the log message `msg` to the detailed log file.
"""
function log_info(msg::String)
    open(detailed_log_file, "a") do f
        println(f, msg)
    end
end

# -------------------------------
# File Naming & Resume Helpers
# -------------------------------

"""
    generate_filename(prefix::String, rec_level::Int)

Generates a filename embedding key parameters (including recursion level and level_initial)
to uniquely identify saved data files.
"""
function generate_filename(prefix::String, rec_level::Int)
    return joinpath(@__DIR__, "data", 
        "$(prefix)_recursion_$(rec_level)_lambda_$(λ_initial)_$(λ_final)_h_$(h_val)_n_$(n_value)_level_$(level_initial)_beta_$(beta_range[1])_$(beta_range[2])_t_$(t_range[1])_$(t_range[2]).txt")
end

"""
    generate_rectangles_filename(rec_level::Int)

Returns the filename for storing rectangle data for the given recursion level.
"""
generate_rectangles_filename(rec_level::Int) = generate_filename("rectangles_with_zeros", rec_level)

"""
    find_saved_rectangles(recursion_steps::Int)

Searches for a saved rectangles file in the data folder matching the current parameters,
with a recursion level lower than `recursion_steps`.
Returns a tuple `(file_path, saved_depth)` if found, or `nothing` otherwise.
"""
function find_saved_rectangles(recursion_steps::Int)
    data_dir = joinpath(@__DIR__, "data")
    files = readdir(data_dir)
    # Allow an optional ".txt" extension
    pattern = Regex("^rectangles_with_zeros_recursion_(\\d+)_lambda_" *
        replace(string(λ_initial), "." => "\\.") * "_" *
        replace(string(λ_final),   "." => "\\.") * "_h_" *
        replace(string(h_val),     "." => "\\.") * "_n_" *
        string(n_value) * "_level_" *
        string(level_initial) * "_beta_" *
        replace(string(beta_range[1]), "." => "\\.") * "_" *
        replace(string(beta_range[2]), "." => "\\.") * "_t_" *
        replace(string(t_range[1]), "." => "\\.") * "_" *
        replace(string(t_range[2]), "." => "\\.") * "(?:\\.txt)?\$")
    best_depth = -1
    best_file = ""
    for f in files
        m = match(pattern, f)
        if m !== nothing
            d = parse(Int, m.captures[1])
            if d > best_depth && d < recursion_steps
                best_depth = d
                best_file = joinpath(data_dir, f)
            end
        end
    end
    if best_file != ""
        return (best_file, best_depth)
    else
        return nothing
    end
end

# -------------------------------
# Model Setup: Hamiltonians & Overlaps
# -------------------------------
include(joinpath(@__DIR__, "params.jl"))
include(joinpath(@__DIR__, "src", "SpinOperatorLibrary.jl"))
include(joinpath(@__DIR__, "src", "LipkinModel.jl"))
using .SpinOperatorLibrary, .LipkinModel

setprecision(BigFloat, precision)
λ_initial = lambda_initial
λ_final   = lambda_final
h_val     = h_value

mkpath(joinpath(@__DIR__, "data"))

# Construct Hamiltonians and compute eigen data.
H_initial = lipkin(λ_initial, h_val, n_value)
H_final   = lipkin(λ_final, h_val, n_value)
e_vals_initial, e_vecs_initial = eigen_sorted(H_initial)
e_vals_final,   e_vecs_final   = eigen_sorted(H_final)

# Form initial state and compute overlap coefficients.
psiI = (e_vecs_initial[:, level_initial] + e_vecs_initial[:, level_initial+1]) / sqrt(BigFloat(2))
dim = length(e_vals_final)
k = [abs(dot(psiI, e_vecs_final[:, i]))^2 for i in 1:dim]

# -------------------------------
# Loschmidt Amplitude & Winding Number
# -------------------------------
"""
    Z(β, t)

Computes the Loschmidt amplitude at inverse temperature `β` and time `t` using high precision.
"""
function Z(β, t)
    s = BigFloat(0) + 0im
    for i in 1:dim
        s += k[i] * exp(-e_vals_final[i] * (β + im*t))
    end
    return s
end

"""
    Z_complex(z::Complex{BigFloat})

Computes the Loschmidt amplitude for a complex argument `z` by using its real and imaginary parts.
"""
Z_complex(z::Complex{BigFloat}) = Z(real(z), imag(z))

"""
    edge_integral(z_start::Complex{BigFloat}, z_end::Complex{BigFloat})

Computes the contour integral along a straight edge from `z_start` to `z_end` for the logarithmic derivative of Z.
"""
function edge_integral(z_start::Complex{BigFloat}, z_end::Complex{BigFloat})
    d = z_end - z_start
    integrand(t) = begin
        z = z_start + t * d
        Z_val = Z_complex(z)
        Z_val = abs(Z_val) < eps(BigFloat) ? eps(BigFloat) : Z_val
        dZ = FiniteDiff.finite_difference_derivative(Z_complex, z, Val(:central), Complex{BigFloat}(1e-10, 1e-10))
        return (dZ / Z_val) * d
    end
    res, _ = quadgk(integrand, BigFloat("0.0"), BigFloat("1.0"); rtol=integration_tolerance, order=integration_order)
    return res
end

"""
    compute_winding_number_rect(zmin::Complex{BigFloat}, zmax::Complex{BigFloat})

Computes the winding number for the rectangular contour defined by corners `zmin` and `zmax`.
"""
function compute_winding_number_rect(zmin::Complex{BigFloat}, zmax::Complex{BigFloat})
    x_min, x_max = real(zmin), real(zmax)
    y_min, y_max = imag(zmin), imag(zmax)
    c1 = Complex{BigFloat}(x_min, y_min)
    c2 = Complex{BigFloat}(x_max, y_min)
    c3 = Complex{BigFloat}(x_max, y_max)
    c4 = Complex{BigFloat}(x_min, y_max)
    total = edge_integral(c1, c2) + edge_integral(c2, c3) + edge_integral(c3, c4) + edge_integral(c4, c1)
    winding = total / (2π * im)
    return real(winding)
end

# -------------------------------
# Subdivision & Recursion
# -------------------------------
"""
    subdivide_rect(zmin::Complex{BigFloat}, zmax::Complex{BigFloat}; threshold=winding_threshold, depth=0)

Subdivides the rectangle defined by `zmin` and `zmax` into smaller subrectangles, computes the winding number
for each, and returns those subrectangles where the absolute winding number exceeds the threshold.
"""
function subdivide_rect(zmin::Complex{BigFloat}, zmax::Complex{BigFloat}; threshold=winding_threshold, depth=0)
    x_min, x_max = real(zmin), real(zmax)
    y_min, y_max = imag(zmin), imag(zmax)
    
    subd_big = BigFloat(subd)
    x_parts = [x_min + (x_max - x_min) * BigFloat(i)/subd_big for i in 0:subd]
    y_parts = [y_min + (y_max - y_min) * BigFloat(i)/subd_big for i in 0:subd]
    
    idxs = [(i, j) for i in 1:subd for j in 1:subd]
    thread_results = [Vector{Tuple{Complex{BigFloat}, Complex{BigFloat}}}() for _ in 1:Threads.nthreads()]
    thread_logs = [String[] for _ in 1:Threads.nthreads()]
    
    Threads.@threads for idx in idxs
        @inbounds begin
            i, j = idx
            tid = Threads.threadid()
            sub_min = Complex{BigFloat}(x_parts[i], y_parts[j])
            sub_max = Complex{BigFloat}(x_parts[i+1], y_parts[j+1])
            wn = compute_winding_number_rect(sub_min, sub_max)
            if abs(wn) > threshold
                push!(thread_logs[tid], "Depth $depth: Accepted [$sub_min, $sub_max] with winding $wn")
                push!(thread_results[tid], (sub_min, sub_max))
            end
        end
    end
    
    for log_buf in thread_logs
        for msg in log_buf
            log_info(msg)
        end
    end
    log_info("Depth $depth: Selected $(sum(length.(thread_results))) of $(subd*subd) subrectangles")
    return vcat(thread_results...)
end

"""
    find_rectangles(zmin::Complex{BigFloat}, zmax::Complex{BigFloat}, depth::Int)

Recursively subdivides the rectangle defined by `zmin` and `zmax` until the maximum recursion depth is reached.
Returns a vector of subrectangles (each represented as a tuple) where zeros are expected.
"""
function find_rectangles(zmin::Complex{BigFloat}, zmax::Complex{BigFloat}, depth::Int)
    log_info("Depth $depth: Processing rectangle [$(zmin), $(zmax)]")
    if depth >= recursion_steps
        log_info("Depth $depth: Maximum recursion reached.")
        return [(zmin, zmax)]
    else
        rects = subdivide_rect(zmin, zmax; depth=depth)
        log_info("Depth $depth: Found $(length(rects)) active subrectangles")
        result = Vector{Tuple{Complex{BigFloat},Complex{BigFloat}}}()
        for (sub_min, sub_max) in rects
            append!(result, find_rectangles(sub_min, sub_max, depth+1))
        end
        log_info("Depth $depth: Returning $(length(result)) rectangles")
        return result
    end
end

"""
    roots_from_rectangles(rectangles)

Computes approximate zeros by taking the midpoint of each rectangle provided.
"""
roots_from_rectangles(rectangles) = [(r[1] + r[2]) / BigFloat(2) for r in rectangles]

# -------------------------------
# Main Data Generation & Resume Logic
# -------------------------------
const zRange = (Complex{BigFloat}(beta_range[1], t_range[1]),
                Complex{BigFloat}(beta_range[2], t_range[2]))
println("Starting search over $zRange ...")
log_info("Starting rectangle search over $zRange")

# Resume from saved rectangles if available.
saved_info = find_saved_rectangles(recursion_steps)
if saved_info !== nothing
    saved_file, saved_depth = saved_info
    println("Found saved file with matching parameters: $saved_file")
    println("Resuming from saved rectangles at recursion depth $saved_depth.")
    log_info("Found saved file $saved_file (depth $saved_depth). Resuming...")
    rects_data = readdlm(saved_file)
    saved_rectangles = [(Complex{BigFloat}(row[1], row[2]), Complex{BigFloat}(row[3], row[4]))
                        for row in eachrow(rects_data)]
    rectangles_with_zeros = saved_depth < recursion_steps ?
        reduce(vcat, [find_rectangles(zmin, zmax, saved_depth) for (zmin, zmax) in saved_rectangles]) :
        saved_rectangles
else
    rectangles_with_zeros = find_rectangles(zRange[1], zRange[2], 0)
end

# Final calculation of the zeros coordinates from rectangles.
approx_roots = [(r[1] + r[2]) / BigFloat(2) for r in rectangles_with_zeros]
println("Total approximate zeros found: $(length(approx_roots))")
log_info("Final approximate roots:")
for root in approx_roots
    log_info("Zero at: $root")
end

elapsed_time = time() - start_time

# Save results with filenames including the recursion level.
rectangles_filename = generate_rectangles_filename(recursion_steps)

writedlm(rectangles_filename, [(real(r[1]), imag(r[1]), real(r[2]), imag(r[2])) for r in rectangles_with_zeros])
writedlm("data/approx_roots.txt", [(real(root), imag(root)) for root in approx_roots])
println("Computation took $(elapsed_time) seconds")