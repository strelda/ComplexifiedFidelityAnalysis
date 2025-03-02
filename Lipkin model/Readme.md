# Lipkin Zeros Finder

Mathematica code served as a proof of concept, which was then rewritten into Julia for efficiency.

Code for finding the zeros of the Loschmidt amplitude in the Lipkin model. It constructs the Lipkin Hamiltonian using spin operator matrices, computes sorted eigenvalues and eigenvectors, and then locates zeros in the complex plane via a contour-integration (winding number) method.

## Project Structure

- **ComplexifiedZeroes.nb:** Mathematica notebook with the original code.
- **library/SpinOperatorLibrary.m:** Mathematica library to construct the spin operator matrices (Jx, Jy, Jz).
- **Project.toml:** Project dependencies.
- **src/SpinOperatorLibrary.jl:** Library to construct the spin operator matrices (Jx, Jy, Jz).
- **src/LipkinModel.jl:** Module to build the Lipkin Hamiltonian and compute sorted eigen systems.
- **params.jl:** File with parameters (model parameters, integration settings, etc.) that can be easily adjusted.

- **find_zeroes.jl:** Main driver script that computes the fidelity function, performs the contour integration to locate zeros, and saves their coordinates.
- **generate_Z_background.jl:** Script to generate the background of the Loschmidt amplitude.
- **plot.jl:** Script to plot the results.

## Running the Project

1. Open a terminal in the project directory.
2. Activate the project environment:
   ```bash
   julia --project=.
   import Pkg
   Pkg.instantiate()
   julia main.jl
    ```

3. Run three scripts to plot the results:
   ```bash
   julia find_zeroes.jl
   julia generate_Z_background.jl
   julia plot.jl
   ```
or run the following command to execute all the scripts at once:
   ```bash
   ./runall.sh
   ```