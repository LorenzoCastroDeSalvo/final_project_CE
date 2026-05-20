using Plots
using HANKPolicies

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "output")

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

function make_policy_value_plots(; output_dir=OUTPUT_DIR)
    mkpath(output_dir)
    p = baseline_params()
    numerical_checks(p)

    a_grid = collect(range(-2.5, 2.5, length=400))
    xi_values = [p.xi_bar - 2.0 * p.sigma, p.xi_bar, p.xi_bar + 2.0 * p.sigma]

    default(fontfamily="Computer Modern", legend=:best)

    plot_against_assets(a_grid, xi_values, c_policy, p,
        joinpath(output_dir, "policy_consumption_vs_assets.png"); ylabel="consumption c")
    plot_against_assets(a_grid, xi_values, labor_policy, p,
        joinpath(output_dir, "policy_labor_vs_assets.png"); ylabel="labor l")
    plot_against_assets(a_grid, xi_values, a_prime_policy, p,
        joinpath(output_dir, "policy_savings_vs_assets.png"); ylabel="next assets a'")

    plot_against_cash_on_hand(a_grid, xi_values, c_policy, p,
        joinpath(output_dir, "policy_consumption_vs_cash_on_hand.png"); ylabel="consumption c")
    plot_against_cash_on_hand(a_grid, xi_values, labor_policy, p,
        joinpath(output_dir, "policy_labor_vs_cash_on_hand.png"); ylabel="labor l")

    plot_against_assets(a_grid, xi_values,
        (a, xi, p) -> flow_utility(c_policy(a, xi, p), labor_policy(a, xi, p), xi, p),
        p, joinpath(output_dir, "flow_utility_vs_assets.png"); ylabel="flow utility")

    value_result = evaluate_value_function(p; n_nodes=15, verbose=true)
    value_xi_values = xi_values
    value_indices = [closest_index(value_result.xi_nodes, xi) for xi in value_xi_values]
    labels = ["xi low", "xi mean", "xi high"]

    plt_value = plot(xlabel="assets a", ylabel="value V(a, xi)", grid=true,
        linewidth=2, size=(760, 500))
    for (idx, label) in zip(value_indices, labels)
        plot!(plt_value, value_result.a_grid, value_result.V[:, idx], label=label)
    end
    savefig(plt_value, joinpath(output_dir, "value_function_vs_assets.png"))

    

    return (; p, value_result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    make_policy_value_plots(output_dir=OUTPUT_DIR)
end
