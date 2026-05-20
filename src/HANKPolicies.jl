module HANKPolicies

using CSV
using DataFrames
using FastGaussQuadrature
using Plots

export HANKParams,
    baseline_params,
    compute_omega,
    cash_on_hand,
    c_policy,
    labor_policy,
    a_prime_policy,
    flow_utility,
    gh_normal_nodes,
    interp_linear_clamped,
    evaluate_value_function,
    numerical_checks,
    make_omega_comparative_statics_plot,
    make_policy_value_plots,
    run

Base.@kwdef struct HANKParams
    R::Float64 = 1.04
    theta::Float64 = 0.85
    wbar::Float64 = 1.0
    ybar::Float64 = 1.0
    gamma::Float64 = 2.0
    rho::Float64 = 1.0 / 3.0
    phi::Float64 = -5.76
    epsilon::Float64 = 10.0
    kappa::Float64 = 0.1
    betatilde::Float64 = theta / R
    mu::Float64 = (1.0 - betatilde) / (1.0 + gamma * rho * wbar)
    sigma::Float64 = 0.5 / (wbar * (1.0 - gamma * rho * mu * wbar))
    Lambda::Float64 = gamma^2 * mu^2 * wbar^2 * sigma^2
    beta::Float64 = exp(-Lambda / 2.0) / R
    Theta::Float64 = 1.0 - Lambda * phi / gamma
    Omega::Float64 = (Lambda + Theta - 1.0) / ((1.0 - betatilde) * (1.0 - Lambda))
    q::Float64 = theta / R
    xi_bar::Float64 = ybar * (1.0 + gamma * rho) - rho * log(wbar)
end

function baseline_params(; R=1.04, theta=0.85, wbar=1.0, ybar=1.0,
    gamma=2.0, rho=1.0 / 3.0, phi=-5.76, epsilon=10.0, kappa=0.1)

    betatilde = theta / R
    mu = (1.0 - betatilde) / (1.0 + gamma * rho * wbar)
    sigma = 0.5 / (wbar * (1.0 - gamma * rho * mu * wbar))
    Lambda = gamma^2 * mu^2 * wbar^2 * sigma^2
    beta = exp(-Lambda / 2.0) / R
    Theta = 1.0 - Lambda * phi / gamma
    Omega = (Lambda + Theta - 1.0) / ((1.0 - betatilde) * (1.0 - Lambda))
    q = theta / R
    xi_bar = ybar * (1.0 + gamma * rho) - rho * log(wbar)

    return HANKParams(; R, theta, wbar, ybar, gamma, rho, phi, epsilon, kappa,
        betatilde, mu, sigma, Lambda, beta, Theta, Omega, q, xi_bar)
end

function compute_omega(; sigma, gamma, rho, phi, R, theta, wbar)
    betatilde = theta / R
    mu = (1.0 - betatilde) / (1.0 + gamma * rho * wbar)
    Lambda = gamma^2 * mu^2 * wbar^2 * sigma^2
    Theta = 1.0 - Lambda * phi / gamma
    Omega = (Lambda + Theta - 1.0) / ((1.0 - betatilde) * (1.0 - Lambda))
    return (; Omega, Lambda, Theta, mu, betatilde)
end

cash_on_hand(a, xi, p::HANKParams) = a + p.wbar * (xi - p.xi_bar)

c_policy(a, xi, p::HANKParams) = p.ybar + p.mu * cash_on_hand(a, xi, p)

function labor_policy(a, xi, p::HANKParams)
    return p.rho * log(p.wbar) - p.gamma * p.rho * c_policy(a, xi, p) + xi
end

function a_prime_policy(a, xi, p::HANKParams)
    c = c_policy(a, xi, p)
    labor = labor_policy(a, xi, p)
    return (p.wbar * labor + a - c) / p.q
end

function flow_utility(c, labor, xi, p::HANKParams)
    return -(1.0 / p.gamma) * exp(-p.gamma * c) - p.rho * exp((labor - xi) / p.rho)
end

function gh_normal_nodes(p::HANKParams; n_nodes=15)
    nodes, weights = gausshermite(n_nodes)
    xi_nodes = p.xi_bar .+ sqrt(2.0) * p.sigma .* nodes
    probs = weights ./ sqrt(pi)
    return xi_nodes, probs
end

function interp_linear_clamped(xs::AbstractVector, ys::AbstractVector, x::Real)
    x <= first(xs) && return first(ys)
    x >= last(xs) && return last(ys)

    idx = searchsortedlast(xs, x)
    x0 = xs[idx]
    x1 = xs[idx + 1]
    y0 = ys[idx]
    y1 = ys[idx + 1]
    weight = (x - x0) / (x1 - x0)
    return (1.0 - weight) * y0 + weight * y1
end

function evaluate_value_function(p::HANKParams;
    a_grid=collect(range(-5.0, 5.0, length=800)),
    n_nodes=15,
    max_iter=10_000,
    tol=1e-8,
    verbose=false)

    xi_nodes, probs = gh_normal_nodes(p; n_nodes)
    n_a = length(a_grid)
    n_xi = length(xi_nodes)
    V = zeros(n_a, n_xi)
    V_new = similar(V)
    discount = p.beta * p.theta

    for iter in 1:max_iter
        for (k, xi) in pairs(xi_nodes)
            for (i, a) in pairs(a_grid)
                c = c_policy(a, xi, p)
                labor = labor_policy(a, xi, p)
                flow_u = flow_utility(c, labor, xi, p)
                a_next = a_prime_policy(a, xi, p)
                continuation = 0.0

                for j in 1:n_xi
                    continuation += probs[j] * interp_linear_clamped(a_grid, view(V, :, j), a_next)
                end

                V_new[i, k] = flow_u + discount * continuation
            end
        end

        diff = maximum(abs.(V_new .- V))
        V, V_new = V_new, V

        if verbose && (iter == 1 || iter % 250 == 0 || diff < tol)
            @info "Policy evaluation" iter diff
        end

        if diff < tol
            return (; V, a_grid, xi_nodes, probs, iterations=iter, supnorm=diff)
        end
    end

    return (; V, a_grid, xi_nodes, probs, iterations=max_iter,
        supnorm=maximum(abs.(V_new .- V)))
end

function numerical_checks(p::HANKParams; atol=1e-8)
    a0, a1 = -0.7, 1.3
    xi0 = p.xi_bar + 0.4 * p.sigma
    x0 = cash_on_hand(a0, xi0, p)
    x1 = cash_on_hand(a1, xi0, p)
    c_slope = (c_policy(a1, xi0, p) - c_policy(a0, xi0, p)) / (x1 - x0)
    labor_slope = (labor_policy(a1, xi0, p) - labor_policy(a0, xi0, p)) / (a1 - a0)

    @assert isapprox(c_slope, p.mu; atol) "c_policy is not linear in cash-on-hand with slope mu"
    @assert isapprox(labor_slope, -p.gamma * p.rho * p.mu; atol) "labor_policy slope in assets is wrong"
    @assert isapprox(a_prime_policy(0.0, p.xi_bar, p), 0.0; atol=1e-10) "a_prime_policy(0, xi_bar) is not zero"
    @assert p.beta * p.theta < 1.0 "beta * theta must be below one"
    @assert p.Lambda < 1.0 "Lambda must be below one for finite Omega"

    return true
end

function closest_index(xs, x)
    return argmin(abs.(xs .- x))
end

function plot_against_assets(a_grid, xi_values, yfun, p, path; ylabel)
    labels = ["xi low", "xi mean", "xi high"]
    plt = plot(xlabel="assets a", ylabel=ylabel, grid=true, linewidth=2, size=(760, 500))
    for (xi, label) in zip(xi_values, labels)
        plot!(plt, a_grid, [yfun(a, xi, p) for a in a_grid], label=label)
    end
    savefig(plt, path)
    return plt
end

function plot_against_cash_on_hand(a_grid, xi_values, yfun, p, path; ylabel)
    labels = ["xi low", "xi mean", "xi high"]
    plt = plot(xlabel="cash-on-hand x", ylabel=ylabel, grid=true, linewidth=2, size=(760, 500))
    for (xi, label) in zip(xi_values, labels)
        x_grid = [cash_on_hand(a, xi, p) for a in a_grid]
        plot!(plt, x_grid, [yfun(a, xi, p) for a in a_grid], label=label)
    end
    savefig(plt, path)
    return plt
end

function omega_sweep_data(p::HANKParams)
    sigmas = collect(range(0.1, 0.9, length=120))
    gammas = collect(range(0.75, 5.0, length=120))
    phis = collect(range(-10.0, -0.25, length=120))

    sigma_df = DataFrame(
        parameter=fill("sigma", length(sigmas)),
        value=sigmas,
        Omega=[compute_omega(; sigma=sigma, gamma=p.gamma, rho=p.rho,
            phi=p.phi, R=p.R, theta=p.theta, wbar=p.wbar).Omega for sigma in sigmas],
    )
    gamma_df = DataFrame(
        parameter=fill("gamma", length(gammas)),
        value=gammas,
        Omega=[compute_omega(; sigma=p.sigma, gamma=gamma, rho=p.rho,
            phi=p.phi, R=p.R, theta=p.theta, wbar=p.wbar).Omega for gamma in gammas],
    )
    phi_df = DataFrame(
        parameter=fill("phi", length(phis)),
        value=phis,
        Omega=[compute_omega(; sigma=p.sigma, gamma=p.gamma, rho=p.rho,
            phi=phi, R=p.R, theta=p.theta, wbar=p.wbar).Omega for phi in phis],
    )

    return vcat(sigma_df, gamma_df, phi_df)
end

function make_omega_comparative_statics_plot(; output_dir="images", p=baseline_params())
    mkpath(output_dir)
    data = omega_sweep_data(p)

    sigma_data = filter(:parameter => ==("sigma"), data)
    gamma_data = filter(:parameter => ==("gamma"), data)
    phi_data = filter(:parameter => ==("phi"), data)

    p_sigma = plot(sigma_data.value, sigma_data.Omega, xlabel="sigma", ylabel="Omega",
        title="Income-risk volatility", legend=false, grid=true, linewidth=2)
    p_gamma = plot(gamma_data.value, gamma_data.Omega, xlabel="gamma", ylabel="Omega",
        title="Risk aversion", legend=false, grid=true, linewidth=2)
    p_phi = plot(phi_data.value, phi_data.Omega, xlabel="phi", ylabel="Omega",
        title="Income-risk cyclicality", legend=false, grid=true, linewidth=2)
    plt = plot(p_sigma, p_gamma, p_phi, layout=(1, 3), size=(1050, 340))

    png_path = joinpath(output_dir, "figure1_omega_comparative_statics.png")
    pdf_path = joinpath(output_dir, "figure1_omega_comparative_statics.pdf")
    csv_path = joinpath(output_dir, "figure1_omega_comparative_statics.csv")
    savefig(plt, png_path)
    savefig(plt, pdf_path)
    CSV.write(csv_path, data)

    return (; plot=plt, data, paths=(png=png_path, pdf=pdf_path, csv=csv_path))
end

function make_policy_value_plots(; output_dir="images",
    value_grid_length=800,
    n_nodes=15,
    value_max_iter=10_000,
    value_tol=1e-8,
    verbose=false)

    mkpath(output_dir)
    p = baseline_params()
    numerical_checks(p)

    a_grid = collect(range(-2.5, 2.5, length=400))
    xi_values = [p.xi_bar - 2.0 * p.sigma, p.xi_bar, p.xi_bar + 2.0 * p.sigma]

    default(fontfamily="Computer Modern", legend=:best)

    paths = Dict{String,String}()
    paths["policy_consumption_vs_assets"] = joinpath(output_dir, "policy_consumption_vs_assets.png")
    plot_against_assets(a_grid, xi_values, c_policy, p,
        paths["policy_consumption_vs_assets"]; ylabel="consumption c")

    paths["policy_labor_vs_assets"] = joinpath(output_dir, "policy_labor_vs_assets.png")
    plot_against_assets(a_grid, xi_values, labor_policy, p,
        paths["policy_labor_vs_assets"]; ylabel="labor l")

    paths["policy_savings_vs_assets"] = joinpath(output_dir, "policy_savings_vs_assets.png")
    plot_against_assets(a_grid, xi_values, a_prime_policy, p,
        paths["policy_savings_vs_assets"]; ylabel="next assets a'")

    paths["policy_consumption_vs_cash_on_hand"] =
        joinpath(output_dir, "policy_consumption_vs_cash_on_hand.png")
    plot_against_cash_on_hand(a_grid, xi_values, c_policy, p,
        paths["policy_consumption_vs_cash_on_hand"]; ylabel="consumption c")

    paths["policy_labor_vs_cash_on_hand"] =
        joinpath(output_dir, "policy_labor_vs_cash_on_hand.png")
    plot_against_cash_on_hand(a_grid, xi_values, labor_policy, p,
        paths["policy_labor_vs_cash_on_hand"]; ylabel="labor l")

    paths["flow_utility_vs_assets"] = joinpath(output_dir, "flow_utility_vs_assets.png")
    plot_against_assets(a_grid, xi_values,
        (a, xi, p) -> flow_utility(c_policy(a, xi, p), labor_policy(a, xi, p), xi, p),
        p, paths["flow_utility_vs_assets"]; ylabel="flow utility")

    value_grid = collect(range(-5.0, 5.0, length=value_grid_length))
    value_result = evaluate_value_function(p; a_grid=value_grid, n_nodes,
        max_iter=value_max_iter, tol=value_tol, verbose)
    value_indices = [closest_index(value_result.xi_nodes, xi) for xi in xi_values]
    labels = ["xi low", "xi mean", "xi high"]

    plt_value = plot(xlabel="assets a", ylabel="value V(a, xi)", grid=true,
        linewidth=2, size=(760, 500))
    for (idx, label) in zip(value_indices, labels)
        plot!(plt_value, value_result.a_grid, value_result.V[:, idx], label=label)
    end
    paths["value_function_vs_assets"] = joinpath(output_dir, "value_function_vs_assets.png")
    savefig(plt_value, paths["value_function_vs_assets"])

    return (; p, value_result, paths)
end

"""
    run(; output_dir="images", kwargs...)

Reproduce the Julia hand-in exhibits for the project and save them in `output_dir`.
The default `run()` is the unique package entry point expected by the CompEcon
project checklist.
"""
function run(; output_dir="images", kwargs...)
    omega = make_omega_comparative_statics_plot(; output_dir)
    policies = make_policy_value_plots(; output_dir, kwargs...)

    return (; output_dir=abspath(output_dir), omega, policies)
end

end
