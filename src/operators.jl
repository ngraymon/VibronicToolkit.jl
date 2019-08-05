# Operators in the harmonic oscillator basis.

"""
Harmonic oscillator basis for `S` surfaces and `M` modes.
"""
struct Basis{S,M}
    "Single-mode basis size."
    size::Int

    "Single-surface basis size."
    dim1::Int
    "Total basis size."
    dim::Int
end

"""
    Basis(::System, size::Int)

Create a basis with `size` basis functions for a single mode.
"""
function Basis(::System{S,M}, size::Int) where {S,M}
    size >= 1 || throw(DomainError(size, "At least 1 basis function."))

    dim1 = size^M
    dim = S * dim1

    Basis{S,M}(size, dim1, dim)
end

"""
    id(basis::Basis)

Identity operator in `basis` for a single mode.
"""
id(basis::Basis) = Matrix{Float64}(I, basis.size, basis.size)

"""
    a(basis::Basis)

Annihilation operator in `basis` for a single mode.
"""
a(basis::Basis) = diagm(1 => sqrt.(1.0:basis.size-1))

"""
    n(basis::Basis)

Number operator in `basis` for a single mode.
"""
n(basis::Basis) = diagm(0 => 0.0:basis.size-1)

"""
    p2(basis::Basis)

Position operator, raised to the power 2, in `basis` for a single mode.
"""
p2(basis::Basis) = -0.5 * (a(basis)' - a(basis))^2

"""
    q(basis::Basis)

Position operator in `basis` for a single mode.
"""
q(basis::Basis) = (a(basis)' + a(basis)) / sqrt(2.0)

"""
    q2(basis::Basis)

Position operator, raised to the power 2, in `basis` for a single mode.
"""
q2(basis::Basis) = q(basis) * q(basis)

"""
    mkop(basis::Basis, op::Matrix{Float64}, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes.
"""
function mkop(basis::Basis{S,M}, op::Matrix{Float64}, idx::Int) where {S,M}
    # Nothing to do for a single mode.
    M == 1 && return op

    ops = []
    for m in 1:M
        if m == idx
            push!(ops, op)
        else
            push!(ops, id(basis))
        end
    end

    result = kron(ops...)
    size(result) == (basis.dim1, basis.dim1) || @warn "Bad size in mkop"

    result
end

"""
    mkop(basis::Basis, op, idx::Int)

Create a multi-mode operator in `basis` that is the tensor product of `op` at
mode `idx` and the identity operator at all other modes, where `op` is a
function that generates the operator in the basis.
"""
mkop(basis::Basis, op, idx::Int) = mkop(basis, op(basis), idx)

"""
    mkid(basis::Basis)

Create a multi-mode identity operator in `basis`.
"""
mkid(basis::Basis) = mkop(basis, id(basis), 1)

"""
    operators(basis::Basis{S,M}, sys::System{S,M})

Create `h0` and `V` operators in `basis` for the system described by `sys`.
"""
function operators(basis::Basis{S,M}, sys::System{S,M}) where {S,M}
    # It's simpler to populate the values into higher rank tensors.
    h0s = zeros(basis.dim1, basis.dim1, S, S)
    Vs = zeros(basis.dim1, basis.dim1, S, S)

    for s1 in 1:S
        for s2 in 1:S
            if s1 == s2
                h0s[:, :, s1, s1] .+= sys.energy[s1, s1] * mkid(basis)

                for m in 1:M
                    h0s[:, :, s1, s1] .+= sys.freq[m, s1] * (mkop(basis, n, m) + 0.5 * mkop(basis, id, m))
                end
            else
                Vs[:, :, s2, s1] .+= sys.energy[s2, s1] * mkid(basis)
            end

            for m in 1:M
                if s1 == s2
                    h0s[:, :, s1, s1] .+= sys.lin[m, s1, s1] * mkop(basis, q, m)
                else
                    Vs[:, :, s2, s1] .+= sys.lin[m, s2, s1] * mkop(basis, q, m)
                end
            end

            for m1 in 1:M
                for m2 in 1:M
                    Vs[:, :, s2, s1] .+= 0.5 * sys.quad[m2, m1, s2, s1] * mkop(basis, q, m2) * mkop(basis, q, m1)
                end
            end
        end
    end

    # Flatten into matrices.
    h0 = reshape(permutedims(h0s, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(h0' - h0)) < 1e-13 || @warn "Asymmetric h0: $(maximum(abs.(h0' - h0)))"
    V = reshape(permutedims(Vs, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(V' - V)) < 1e-13 || @warn "Asymmetric V: $(maximum(abs.(V' - V)))"

    h0, V
end

"""
    operators_h_U(basis::Basis{S,M}, sys::System{S,M})

Create `h0` and `U` operators in `basis` for the system described by `sys`.
Use the explicity operator forms, which is less accurate than the harmonic (n+05) form but is necessary to compare to the T+V splitting.
"""
function operators_h_U(basis::Basis{S,M}, sys::System{S,M}) where {S,M}
    # It's simpler to populate the values into higher rank tensors.
    h0s = zeros(basis.dim1, basis.dim1, S, S)
    Us = zeros(basis.dim1, basis.dim1, S, S)

    for s1 in 1:S
        for s2 in 1:S
            if s1 == s2
                h0s[:, :, s1, s1] .+= sys.energy[s1, s1] * mkid(basis)

                for m in 1:M
                    h0s[:, :, s1, s1] .+= sys.freq[m, s1] * 0.5 * mkop(basis, p2, m)
                    h0s[:, :, s1, s1] .+= sys.freq[m, s1] * 0.5 * mkop(basis, q2, m)
                end
            else
                Us[:, :, s2, s1] .+= sys.energy[s2, s1] * mkid(basis)
            end

            for m in 1:M
                if s1 == s2
                    h0s[:, :, s1, s1] .+= sys.lin[m, s1, s1] * mkop(basis, q, m)
                else
                    Us[:, :, s2, s1] .+= sys.lin[m, s2, s1] * mkop(basis, q, m)
                end
            end

            for m1 in 1:M
                for m2 in 1:M
                    Us[:, :, s2, s1] .+= 0.5 * sys.quad[m2, m1, s2, s1] * mkop(basis, q, m2) * mkop(basis, q, m1)
                end
            end
        end
    end

    # Flatten into matrices.
    h0 = reshape(permutedims(h0s, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(h0' - h0)) < 1e-13 || @warn "Asymmetric h0: $(maximum(abs.(h0' - h0)))"
    U = reshape(permutedims(Us, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(U' - U)) < 1e-13 || @warn "Asymmetric U: $(maximum(abs.(U' - U)))"

    h0, U
end

"""
    operators_T_V(basis::Basis{S,M}, sys::System{S,M})

Create `T` and `V` operators in `basis` for the system described by `sys`.
Use the explicity operator forms, which is less accurate than the harmonic (n+05) form but is necessary to compare to the T+V splitting.
"""
function operators_T_V(basis::Basis{S,M}, sys::System{S,M}) where {S,M}
    # It's simpler to populate the values into higher rank tensors.
    Ts = zeros(basis.dim1, basis.dim1, S, S)
    Vs = zeros(basis.dim1, basis.dim1, S, S)

    for s1 in 1:S
        for s2 in 1:S
            if s1 == s2
                Vs[:, :, s1, s1] .+= sys.energy[s1, s1] * mkid(basis)

                for m in 1:M
                    Ts[:, :, s1, s1] .+= sys.freq[m, s1] * 0.5 * mkop(basis, p2, m)
                    Vs[:, :, s1, s1] .+= sys.freq[m, s1] * 0.5 * mkop(basis, q2, m)
                end
            else
                Vs[:, :, s2, s1] .+= sys.energy[s2, s1] * mkid(basis)
            end

            for m in 1:M
                if s1 == s2
                    Vs[:, :, s1, s1] .+= sys.lin[m, s1, s1] * mkop(basis, q, m)
                else
                    Vs[:, :, s2, s1] .+= sys.lin[m, s2, s1] * mkop(basis, q, m)
                end
            end

            for m1 in 1:M
                for m2 in 1:M
                    Vs[:, :, s2, s1] .+= 0.5 * sys.quad[m2, m1, s2, s1] * mkop(basis, q, m2) * mkop(basis, q, m1)
                end
            end
        end
    end

    # Flatten into matrices.
    T = reshape(permutedims(Ts, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(T' - T)) < 1e-13 || @warn "Asymmetric T: $(maximum(abs.(T' - T)))"
    V = reshape(permutedims(Vs, [1, 3, 2, 4]), (basis.dim, basis.dim))
    maximum(abs.(V' - V)) < 1e-13 || @warn "Asymmetric V: $(maximum(abs.(V' - V)))"

    T, V
end

"""
    vectors(basis::Basis{S,M})

Generate an array of basis vectors for `basis` in the order corresponding to
the layout of the operator matrices.

Each column is a basis vector, with the first `M` values giving the (0-indexed)
state labels of the modes, and the final value giving the (1-indexed) surface
label.
"""
function vectors(basis::Basis{S,M}) where {S,M}
    result = Array{Int}(undef, M+1, basis.dim)

    conf = zeros(Int, M)
    idx = 1
    for s in 1:S
        first = true
        while true
            for m in M:-1:1
                !first && (conf[m] += 1)
                if conf[m] < basis.size
                    break
                else
                    conf[m] = 0
                end
            end
            if first
                first = false
            elseif iszero(conf)
                break
            end

            result[1:M, idx] .= conf
            result[M+1, idx] = s

            idx += 1
        end
    end
    idx-1 == basis.dim || @warn "Invalid number of basis vectors: $(idx-1)"

    result
end

"""
    ptrace_modes(basis::Basis, rho::AbstractMatrix{Float64})

Perform the partial trace over the modes for `rho` in `basis`.
"""
ptrace_modes(basis::Basis{S,M}, rho::AbstractMatrix{Float64}) where {S,M} = ptrace(rho, (basis.dim1, S), 1)
