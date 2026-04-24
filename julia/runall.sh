#!/bin/bash

# Activate the Julia project environment
julia --project=. -e 'import Pkg; Pkg.instantiate()'

# Run the Julia scripts
julia --project=. -t auto find_zeroes.jl
julia --project=. -t auto generate_Z_background.jl
julia --project=. -t auto plot.jl

echo "All scripts executed."
