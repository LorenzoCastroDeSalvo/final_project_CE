using Test
using HANKPolicies

@testset "HANKPolicies package" begin
    p = baseline_params()

    @test numerical_checks(p)
    @test isfinite(p.Omega)
    @test p.beta * p.theta < 1
    @test c_policy(0.0, p.xi_bar, p) ≈ p.ybar
    @test a_prime_policy(0.0, p.xi_bar, p) ≈ 0.0 atol = 1e-10

    out = mktempdir()
    result = HANKPolicies.run(;
        output_dir=out,
        value_grid_length=41,
        n_nodes=5,
        value_max_iter=2,
        value_tol=Inf,
    )

    @test result.output_dir == abspath(out)
    @test isfile(joinpath(out, "figure1_omega_comparative_statics.png"))
    @test isfile(joinpath(out, "figure1_omega_comparative_statics.csv"))
    @test isfile(joinpath(out, "policy_consumption_vs_assets.png"))
    @test isfile(joinpath(out, "value_function_vs_assets.png"))
end
