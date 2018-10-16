# Trotter factorized solution for systems with few modes.

"""
Trotter factorized solution for a small system.
"""
abstract type AbstractTrotter <: Solution end

"""
Trotter factorized solution for a small system at finite temperature.
"""
struct Trotter <: AbstractTrotter
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    Trotter(sys::System, beta::Float64, basis_size::Int, P::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions
and `P` links.
"""
function Trotter(sys::System, beta::Float64, basis_size::Int, P::Int)
    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    # Single Trotter product.
    tau = beta / P
    rho = exp(-tau * h0) * exp(-tau * V)

    # Full path.
    path = rho^P

    Z = tr(path)
    E = tr(path * (h0 + V)) / Z
    Cv = (tr(path * (h0 + V)^2) / Z - E^2) * beta^2

    Trotter(Z, E, Cv)
end

"""
Trotter factorized solution for a small PIGS system.
"""
struct PigsTrotter <: AbstractTrotter
    "Pseudo-partition function."
    Z::Float64
    "Energy."
    E::Float64
end

"""
    PigsTrotter(sys::System{S,M}, trial::TrialWavefunction{S,M}, beta::Float64, basis_size::Int, P::Int)

Calculate the solution for `sys` with `trial` propagated by `beta` using
`basis_size` basis functions and `P` links.
"""
function PigsTrotter(sys::System{S,M}, trial::TrialWavefunction{S,M}, beta::Float64, basis_size::Int, P::Int) where {S,M}
    P % 2 == 0 || throw(DomainError(P, "Number of links must be even."))

    basis = Basis(sys, basis_size)
    h0, V = operators(basis, sys)

    # Single Trotter product.
    tau = beta / P
    rho = exp(-0.5 * tau * V) * exp(-tau * h0) * exp(-0.5 * tau * V)

    # Full path.
    path = rho^P

    trial_vec = trial_mode(trial, basis)

    Z = dot(trial_vec, path * trial_vec)
    E = dot(trial_vec, path * (h0 + V) * trial_vec) / Z

    PigsTrotter(Z, E)
end
