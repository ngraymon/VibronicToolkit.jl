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
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]
basis_size = c[:basis_size]

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    sos = PigsSumOverStates(sys, trial, beta, basis_size)

    # these 'sos' results represent a finite beta result
    println("Z_sos: $(sos.Z)")
    println("E_sos: $(sos.E)")
    println("SvN_sos: $(sos.SvN)")
    println("S2_sos: $(sos.S2)")

    # these 'exact' represent the infinite beta results (i.e. T = 0 Kelvin)
    println("E0_exact: $(sos.E0_exact)")
    println("SvN_exact: $(sos.SvN_exact)")
    println("S2_exact: $(sos.S2_exact)")

	println("beta: $(beta)")
else
    sos = SumOverStates(sys, beta, basis_size)

    println("Z_sos: $(sos.Z)")
    println("E_sos: $(sos.E)")
    println("Cv_sos: $(sos.Cv)")
	println("beta: $(beta)")
end
