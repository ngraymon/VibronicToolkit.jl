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
    "--uncoupled"
        help = "remove inter-surface and quadratic coupling from system"
        action = :store_true
end
c = parse_args(ARGS, s, as_symbols=true)

sys = read(c[:conf], DenseSystem)
beta = c[:beta]

if c[:uncoupled]
    new_sys = simplify(diag(sys))
else
    try
        global new_sys = DiagonalSystem(sys)
    catch ex
        ex isa SurfaceCouplingException || rethrow()
        error("Analytical solution only applies to uncoupled systems")
    end
end

delta_beta = 2.0E-4
beta_plus = beta + 2.0E-4
beta_minus = beta - 2.0E-4

analytical = Analytical(simplify(sys), beta)
analytical_plus = Analytical(simplify(sys), beta_plus)
analytical_minus = Analytical(simplify(sys), beta_minus)

# println("Z: $(analytical.Z)")
# println("E: $(analytical.E)")
# println("Cv: $(analytical.Cv)")
# println("ZH: $(analytical.Z)")
println("Zrho: $(analytical.Z)")
println("Zrho+(beta): $(analytical.Z)")
println("Zrho-(beta): $(analytical.Z)")
# println("E: $(analytical.E)")
println("Erho: $(analytical.E)")
# println("Cv: $(analytical.Cv)")
println("Cvrho: $(analytical.Cv)")
println("beta: $(beta)")

