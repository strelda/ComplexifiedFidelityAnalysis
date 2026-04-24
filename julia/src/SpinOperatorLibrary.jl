module SpinOperatorLibrary

using LinearAlgebra, SparseArrays

export je, Jx, Jy, Jz, Jx2

"""
    je(sn, j, m)

Calculates the matrix element for spin operators.
"""
function je(sn, j, m)
    return sqrt(j*(j + 1) - m*(m + sn))
end

"""
    Jx(n::Int)

Generates the Jx matrix for a spin of j = n/2 (resulting in an (n+1)×(n+1) matrix).
"""
function Jx(n::Int)
    j = n / 2
    N = n + 1
    row = Int[]
    col = Int[]
    vals = Float64[]
    # Upper diagonal: for i = 1 to n, with m = j - i
    for i in 1:n
        m = j - i
        push!(row, i)
        push!(col, i + 1)
        push!(vals, je(1, j, m) / 2)
    end
    # Lower diagonal: same as above (real numbers, so conjugate is the same)
    for i in 1:n
        m = j - i
        push!(row, i + 1)
        push!(col, i)
        push!(vals, je(1, j, m) / 2)
    end
    return sparse(row, col, vals, N, N)
end

"""
    Jy(n::Int)

Generates the Jy matrix for a spin of j = n/2 (resulting in an (n+1)×(n+1) matrix).
"""
function Jy(n::Int)
    j = n / 2
    N = n + 1
    row = Int[]
    col = Int[]
    vals = ComplexF64[]
    # Upper diagonal: using factor 1/(2im)
    for i in 1:n
        m = j - i
        push!(row, i)
        push!(col, i + 1)
        push!(vals, je(1, j, m) / (2im))
    end
    # Lower diagonal: note that the complex conjugate of 1/(2im) is -1/(2im)
    for i in 1:n
        m = j - i
        push!(row, i + 1)
        push!(col, i)
        push!(vals, -je(1, j, m) / (2im))
    end
    return sparse(row, col, vals, N, N)
end

"""
    Jz(n::Int)

Generates the Jz matrix (diagonal) for a spin of j = n/2 (resulting in an (n+1)×(n+1) matrix).
"""
function Jz(n::Int)
    j = n / 2
    N = n + 1
    diag = [j - (i - 1) for i in 1:N]
    return spdiagm(0 => diag)
end

end # module
