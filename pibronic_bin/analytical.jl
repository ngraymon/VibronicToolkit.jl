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
    "--uncoupled"
        help = "remove inter-surface and quadratic coupling from system"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]

if c[:uncoupled]
    sys = simplify(diag(sys))
else
    try
        global sys = DiagonalSystem(sys)
    catch ex
        ex isa SurfaceCouplingException || rethrow()
        error("Analytical solution only applies to uncoupled systems")
    end
end

if c[:pigs]
    trial = UniformTrialWavefunction(sys)
    analytical = PigsAnalytical(sys, trial, beta)

    println("Z_analytical: $(analytical.Z)")
    println("E_analytical: $(analytical.E)")
    println("SvN_analytical: $(analytical.SvN)")
    println("S2_analytical: $(analytical.S2)")
	println("beta: $(beta)")
else
    analytical = Analytical(sys, beta)

    println("Z_analytical: $(analytical.Z)")
    println("E_analytical: $(analytical.E)")
    println("Cv_analytical: $(analytical.Cv)")
	println("beta: $(beta)")
end
