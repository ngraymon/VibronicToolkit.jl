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
    Trotter_h_U(sys::System, beta::Float64, basis_size::Int, P::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions
and `P` links.
Special version which uses the explicit operator forms instead of the more accurate harmonic oscillator (n + 1/2) version.
This uses the same splitting of the Hamiltonian as in SumOverStates, into the harmonic part (h), and everything else (U).
This is for comparison with the SumOverStates_T+V function.
"""
function Trotter_h_U(sys::System, beta::Float64, basis_size::Int, P::Int)
    basis = Basis(sys, basis_size)
    h0, U = operators_h_U(basis, sys)

    # Single Trotter product.
    tau = beta / P
    rho = exp(-tau * h0) * exp(-tau * U)

    # Full path.
    path = rho^P

    Z = tr(path)
    E = tr(path * (h0 + U)) / Z
    Cv = (tr(path * (h0 + U)^2) / Z - E^2) * beta^2

    Trotter(Z, E, Cv)
end

"""
    Trotter_T+V(sys::System, beta::Float64, basis_size::Int, P::Int)

Calculate the solution for `sys` at `beta` with `basis_size` basis functions
and `P` links.
Special version which uses the explicit operator forms.
This uses a traditional splitting of the Hamiltonian into the kinetic energy (T) and everything else (V).
This is for comparison with the SumOverStates_h+U function.
"""
function Trotter_T_V(sys::System, beta::Float64, basis_size::Int, P::Int)
    basis = Basis(sys, basis_size)
    T, V = operators_T_V(basis, sys)

    # Single Trotter product.
    tau = beta / P
    rho = exp(-tau * T) * exp(-tau * V)

    # Full path.
    path = rho^P

    Z = tr(path)
    E = tr(path * (T + V)) / Z
    Cv = (tr(path * (T + V)^2) / Z - E^2) * beta^2

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
    "Von Neumann entanglement entropy."
    SvN::Float64
    "Order-2 Rényi entanglement entropy."
    S2::Float64
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
    prop = exp(-0.5 * tau * V) * exp(-tau * h0) * exp(-0.5 * tau * V)

    # Half path and full path.
    path_half = prop^div(P, 2)
    path = path_half^2

    # Trial wavefunction.
    trial_vec = trial_mode(trial, basis)

    # Density matrix and reduced density matrix.
    rho = (path_half * trial_vec) * (trial_vec' * path_half)
    rho_surface = ptrace_modes(basis, rho)

    Z = dot(trial_vec, path * trial_vec)
    E = dot(trial_vec, path * (h0 + V) * trial_vec) / Z
    SvN = S_vn(rho_surface / Z)
    S2 = S_renyi(rho_surface / Z)

    PigsTrotter(Z, E, SvN, S2)
end
