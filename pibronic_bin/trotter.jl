#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "--pigs"
        help = "use PIGS instead of finite temperature"
        action = :store_true
    "--conf"
        metavar = "FILE"
        help = "path to config file"
        required = true
    "--beta"
        metavar = "T"
        help = "reciprocal temperature"
        arg_type = Float64
        required = true
    "--basis-size"
        metavar = "N"
        help = "single-mode basis size"
        arg_type = Int
        required = true
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
        required = true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
basis_size = c[:basis_size]
P = c[:num_links]

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    trotter = PigsTrotter(sys, trial, beta, basis_size, P)

    println("Z_trotter: $(trotter.Z)")
    println("E_trotter: $(trotter.E)")
    println("SvN_trotter: $(trotter.SvN)")
    println("S2_trotter: $(trotter.S2)")
	println("beta: $(beta)")
	println("tau: $(beta / (P-1))")
else
    trotter = Trotter(sys, beta, basis_size, P)

    println("Z_trotter: $(trotter.Z)")
    println("E_trotter: $(trotter.E)")
    println("Cv_trotter: $(trotter.Cv)")
	println("beta: $(beta)")
	println("tau: $(beta / P)")
end
