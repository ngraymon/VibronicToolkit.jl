let sys = DenseSystem(test_params...),
    beta = 12.34,
    basis_size = 10,
    sos = SumOverStates(sys, beta, basis_size)

    @test isapprox(sos.Z, 7.74822371e-3)
    @test isapprox(sos.E, 3.93864802e-1)
    @test isapprox(sos.Cv, 1.04927608e-10)
end

let sys = DenseSystem(test_params...),
    trial = UniformTrialWavefunction(sys),
    beta = 12.34,
    basis_size = 10,
    sos = PigsSumOverStates(sys, trial, beta, basis_size)

    @test isapprox(sos.Z, 2.69265083e-1)
    @test isapprox(sos.E, 3.93864802e-1)
    @test isapprox(sos.SvN, 9.46869274e-1)
    @test isapprox(sos.S2, 8.27539873e-1)
    @test isapprox(sos.E0_exact, 3.93864802e-1)
    @test isapprox(sos.SvN_exact, 9.46869260e-1)
    @test isapprox(sos.S2_exact, 8.27539848e-1)
end
