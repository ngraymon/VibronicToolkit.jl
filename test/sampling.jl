using VibronicToolkit: SamplingParameters

let sys = DenseSystem(test_params...),
    beta = 12.34,
    P = 256,
    sp = SamplingParameters(simplify(diag(sys)), beta, P)

    @test sp.tau == beta/P
    @test isapprox(weights(sp), [0.999996583, 3.41696734e-6, 1.16757213e-11])
    @test isapprox(weights(sp), weights(simplify(diag(sys)), beta))
end
