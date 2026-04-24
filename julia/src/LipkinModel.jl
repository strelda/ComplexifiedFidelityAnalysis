module LipkinModel

using LinearAlgebra, SparseArrays, ..SpinOperatorLibrary
using GenericLinearAlgebra

export lipkin, eigen_sorted, eigenvalues_sorted, eigenvectors_sorted

"""
Constructs the Lipkin Hamiltonian:
    H = -λ/n * (Jx(n) * Jx(n)) + h * Jz(n)
where n is an even integer (matrix dimension = n + 1).
"""
function lipkin(λ, h, n::Int)
    H = -λ / n * (Jx(n) * Jx(n)) + h * Jz(n)
    # Return the real part (to remove possible tiny imaginary numerical noise)
    return real(H)
end

"""
Computes and returns the eigenvalues and normalized eigenvectors (sorted in ascending order).
Automatically uses a generic eigen solver for BigFloat matrices.
"""
function eigen_sorted(H)
    if eltype(H) <: BigFloat
        # Convert sparse matrix to dense if needed.
        A = issparse(H) ? Matrix(H) : H
        # Use GenericLinearAlgebra for BigFloat matrices.
        F = GenericLinearAlgebra.eigen(A)
        idx = sortperm(F.values)
        sorted_vals = F.values[idx]
        sorted_vecs = F.vectors[:, idx]
        return sorted_vals, sorted_vecs
    else
        A = Matrix(H)  # Ensure H is dense.
        vals, vecs = eigen(A)
        idx = sortperm(vals)
        sorted_vals = vals[idx]
        sorted_vecs = vecs[:, idx]
        return sorted_vals, sorted_vecs
    end
end


eigenvalues_sorted(H) = eigen_sorted(H)[1]
eigenvectors_sorted(H) = eigen_sorted(H)[2]

end # module