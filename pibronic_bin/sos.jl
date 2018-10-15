#!/usr/bin/env julia

using VibronicToolkit

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
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

sos = SumOverStates(sys, beta, basis_size)

println("Z_sos: $(sos.Z)")
println("E_sos: $(sos.E)")
println("Cv_sos: $(sos.Cv)")
println("beta: $(beta)")
