### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# â•”â•â•¡ d2b17e00-382c-11f1-3b3e-63fa5a148559
begin 
	
############################################################
# FIGURE 1 â€” Dogra et al.
# Comparative statics of Î© with respect to Ïƒ, Î³, and Ï•
#
# Required libraries:
#   pkg> add Roots Plots
############################################################

using Roots
using Plots

# plotting backend
gr()

end

# â•”â•â•¡ d3eec6f9-49fb-4b86-a695-a424aaeaffe1
using LinearAlgebra

# â•”â•â•¡ de9a67a3-c904-45eb-a46d-506691dce3a4


# â•”â•â•¡ b55b6eeb-a5df-4367-a64d-420323f32247

# â•”â•â•¡ 9d092bd0-53c6-11f1-2b21-b7063e2b89bd
begin
    PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
    FIGURES_DIR = joinpath(PROJECT_ROOT, "figures")
    mkpath(FIGURES_DIR)

    function save_figure(fig, name::AbstractString)
        png_path = joinpath(FIGURES_DIR, name * ".png")
        pdf_path = joinpath(FIGURES_DIR, name * ".pdf")

        savefig(fig, png_path)
        savefig(fig, pdf_path)

        println("Saved: " * png_path)
        println("Saved: " * pdf_path)
    end
end


# â•”â•â•¡ 81ade795-8250-4957-8e16-ee7bb42242a6
begin

    ############################################################
    # 1. Baseline calibration
    #
    # This cell reproduces Figure 1 by computing Î© on three
    # one-dimensional grids:
    #   (a) varying ÏƒÌ„
    #   (b) varying Î³, while plotting against the implied Î³-grid
    #   (c) varying Ï•
    #
    # The parameter values follow the baseline calibration used
    # in the MATLAB replication code.
    ############################################################

    theta0    = 1.0 - 0.15
    R0        = 1.04
    gambar0   = 2.0
    rhobar0   = 1.0 / 3.0
    sigy0     = 0.5
    dsigydy0  = -3.0

    ############################################################
    # 2. Helper functions
    ############################################################

    # Solve for w from the nonlinear steady-state condition.
    #
    # The function first tries a local root-finding step using the
    # previous solution as a starting value. If that fails, it falls
    # back to an automatic bracketing procedure over a wide grid.
    # This makes the code more robust when parameters move across
    # the comparative-statics grids.
    function solve_w(ec; guess=10.0)
        f(w) = begin
            if !(isfinite(w) && w > 0)
                return NaN
            end

            gr   = ec.gambar * ec.rhobar
            btil = ec.theta / ec.R

            denom = 1.0 + gr - ec.rhobar * log(w)
            if !(isfinite(denom) && denom > 0)
                return NaN
            end

            y   = (1.0 + gr) / denom
            mu  = (1.0 - btil) / (1.0 + gr * w)
            gam = ec.gambar / y
            sig = ec.sigbar * exp(ec.varphi * (y - 1.0))

            LAM = (gam * mu * w * sig)^2
            if !(isfinite(LAM) && abs(1.0 - LAM) > 1e-12)
                return NaN
            end

            THETA = 1.0 - LAM * ec.varphi / gam
            OMEGA_code = (THETA - 1.0 + LAM) / (1.0 - LAM)

            # Nonlinear equilibrium condition corresponding to optw.m
            OMEGA_code - (1.0 - btil) * (w - 1.0) / (1.0 + gr * w)
        end

        # First attempt: fast local solve near the previous root.
        try
            w = find_zero(f, guess, Order1())
            if isfinite(w) && w > 0
                return w
            end
        catch
        end

        # Fallback: search for a sign change on a wide exponential grid.
        grid = exp.(range(log(0.02), log(50.0), length=1200))
        vals = [f(w) for w in grid]

        for i in 1:length(grid)-1
            fa, fb = vals[i], vals[i+1]
            if isfinite(fa) && isfinite(fb) && sign(fa) != sign(fb)
                return find_zero(f, (grid[i], grid[i+1]), Bisection())
            end
        end

        error("No root found for w.")
    end

    # Construct the baseline economic environment and the derived
    # steady-state objects needed for the comparative statics.
    function makeec_baseline(; theta=theta0, R=R0, gambar=gambar0, rhobar=rhobar0,
                               sigy=sigy0, dsigydy=dsigydy0)

        gr    = gambar * rhobar
        btil  = theta / R
        sigbar = (1.0 + gr) * sigy / (1.0 + btil * gr)

        if sigy > 0
            varphi = dsigydy / sigy + gambar * (1.0 - btil) / (1.0 + btil * gr)
        else
            varphi = 0.0
        end

        ec = (
            theta  = theta,
            R      = R,
            gambar = gambar,
            rhobar = rhobar,
            sigy   = sigy,
            dsigydy = dsigydy,
            gr     = gr,
            btil   = btil,
            sigbar = sigbar,
            varphi = varphi
        )

        w = solve_w(ec)

        merge(ec, (w=w,))
    end

    # Given a parameterized environment, compute the steady-state
    # statistics used in the figure.
    #
    # OMEGA_code is the object naturally produced by the internal
    # MATLAB formulas, while Omega rescales it into the Î© object
    # shown in the paper.
    function make_compstats_like(ec; guess=10.0)
        w = solve_w(ec; guess=guess)

        gr   = ec.gambar * ec.rhobar
        btil = ec.theta / ec.R

        y   = (1.0 + gr) / (1.0 + gr - ec.rhobar * log(w))
        mu  = (1.0 - btil) / (1.0 + gr * w)
        gam = ec.gambar / y
        sig = ec.sigbar * exp(ec.varphi * (y - 1.0))

        LAM = (gam * mu * w * sig)^2
        THETA = 1.0 - LAM * ec.varphi / gam

        # Internal object from the replication-package formulas.
        OMEGA_code = (THETA - 1.0 + LAM) / (1.0 - LAM)

        # Î© plotted in Figure 1 of the paper.
        Omega = OMEGA_code / (1.0 - btil)

        (
            w = w,
            y = y,
            mu = mu,
            gam = gam,
            sig = sig,
            LAM = LAM,
            THETA = THETA,
            OMEGA = OMEGA_code,
            Omega = Omega
        )
    end

    ############################################################
    # 3. Baseline object
    ############################################################

    ec0 = makeec_baseline()

    ############################################################
    # 4. Grids
    #
    # These grids follow the MATLAB replication logic:
    #   - ÏƒÌ„ grid for panel (a)
    #   - Ï• grid for panel (c)
    #   - Î³Ì„ input grid for panel (b), later mapped into gam_grid
    ############################################################

    n = 50
    sig_grid    = collect(range(0.001 * ec0.sigbar, 1.0, length=n))
    varphi_grid = collect(range(-6.0, 6.0, length=n))
    gambar_grid = collect(range(0.001, 5.0, length=n))

    ############################################################
    # 5. Figure 1 objects
    ############################################################

    Omega_sig_grid    = zeros(n)
    Omega_varphi_grid = zeros(n)
    Omega_gambar_grid = zeros(n)
    gam_grid          = zeros(n)

    # Panel (a): vary ÏƒÌ„.
    # The previous solution is reused as the next initial guess
    # because neighboring grid points have very similar roots.
    local wguess = ec0.w
    for i in eachindex(sig_grid)
        ecuse = merge(ec0, (sigbar = sig_grid[i],))
        out = make_compstats_like(ecuse; guess=wguess)
        Omega_sig_grid[i] = out.Omega
        wguess = out.w
    end

    # Panel (c): vary Ï•.
    wguess = ec0.w
    for i in eachindex(varphi_grid)
        ecuse = merge(ec0, (varphi = varphi_grid[i],))
        out = make_compstats_like(ecuse; guess=wguess)
        Omega_varphi_grid[i] = out.Omega
        wguess = out.w
    end

    # Panel (b): vary Î³Ì„ in the primitive parametrization, but
    # plot the result against the implied Î³ values returned by
    # the steady-state computation.
    wguess = ec0.w
    for i in eachindex(gambar_grid)
        ecuse = merge(ec0, (gambar = gambar_grid[i],))
        out = make_compstats_like(ecuse; guess=wguess)
        Omega_gambar_grid[i] = out.Omega
        gam_grid[i] = out.gam
        wguess = out.w
    end

    ############################################################
    # 6. Plot
    ############################################################

    default(
        legend=false,
        grid=true,
        gridalpha=0.30,
        gridlinewidth=0.7,
        gridcolor=:lightgray,
        foreground_color_axis=:black,
        foreground_color_border=:black,
        background_color=:white,
        framestyle=:box,
        linewidth=2.0,
        size=(1500, 400),
        titlefont=font(20, "Computer Modern"),
        guidefont=font(18, "Computer Modern"),
        tickfont=font(11, "Computer Modern")
    )

    p1 = plot(
        sig_grid, Omega_sig_grid,
        color=:blue,
        title="(a)",
        xlabel="Ïƒ",
        ylabel="Î©",
        xlims=(sig_grid[1], sig_grid[end]),
        ylims=(0.0, 0.40),
        yticks=0.0:0.05:0.40
    )

    p2 = plot(
        gam_grid, Omega_gambar_grid,
        color=:blue,
        title="(b)",
        xlabel="Î³",
        ylabel="",
        xlims=(gam_grid[1], gam_grid[end]),
        ylims=(0.0, 0.25),
        yticks=0.0:0.05:0.25
    )

    p3 = plot(
        varphi_grid, Omega_varphi_grid,
        color=:blue,
        title="(c)",
        xlabel="Ï•",
        ylabel="",
        xlims=(varphi_grid[1], varphi_grid[end]),
        ylims=(-0.10, 0.20),
        yticks=-0.10:0.05:0.20
    )

    hline!(p3, [0.0], color=:black, lw=1.0)
    fig1 = plot(p1, p2, p3, layout=(1, 3), margin=5Plots.mm)

    save_figure(fig1, "figure1")

    fig1
end

# â•”â•â•¡ 1f61e41a-3382-474d-b0d5-994933f5b722
begin
    

    ############################################################
    # FIGURE 2 â€” Dogra et al.
    # The effect of monetary policy on consumption inequality Î£_t
    ############################################################

    baseline_fig2 = (
        theta = 1.0 - 0.15,
        R = 1.04,
        gambar = 2.0,
        rhobar = 1.0 / 3.0,
        sigy = 0.5,
        dsigydy = -3.0,
        varphi_override = nothing
    )

    # Build the parameters that are implied by the baseline calibration
    function derived_params_fig2(p)
        gr      = p.gambar * p.rhobar
        btil    = p.theta / p.R
        sigbar  = (1.0 + gr) * p.sigy / (1.0 + btil * gr)

        correction = p.gambar * (1.0 - btil) / (1.0 + btil * gr)
        varphi = isnothing(p.varphi_override) ? p.dsigydy / p.sigy + correction : p.varphi_override

        return (; gr, btil, sigbar, varphi)
    end

    # Nonlinear condition used to pin down w in steady state
    function omega_residual_fig2(w, p)
        d = derived_params_fig2(p)

        y   = (1.0 + d.gr) / (1.0 + d.gr - p.rhobar * log(w))
        gam = p.gambar / y
        mu  = (1.0 - d.btil) / (1.0 + d.gr * w)
        sig = d.sigbar * exp(d.varphi * (y - 1.0))

        LAM   = (gam * mu * w * sig)^2
        THETA = 1.0 - LAM * d.varphi / gam
        OMEGA = (THETA - 1.0 + LAM) / (1.0 - LAM)

        return OMEGA - (1.0 - d.btil) * (w - 1.0) / (1.0 + d.gr * w)
    end

    # Solve the steady state and collect the objects needed for the IRFs
    function steady_state_fig2(p)
        d = derived_params_fig2(p)

        w = find_zero(w -> omega_residual_fig2(w, p), (0.2, 25.0), Bisection())

        y   = (1.0 + d.gr) / (1.0 + d.gr - p.rhobar * log(w))
        gam = p.gambar / y
        rho = p.rhobar * y
        mu  = (1.0 - d.btil) / (1.0 + d.gr * w)
        sig = d.sigbar * exp(d.varphi * (y - 1.0))

        LAM   = (gam * mu * w * sig)^2
        THETA = 1.0 - LAM * d.varphi / gam
        OMEGA = (THETA - 1.0 + LAM) / (1.0 - LAM)

        bet = (1.0 / p.R) * exp(-LAM / 2.0)
        sigma_lag_coeff = d.btil / bet

        return (
            w = w,
            y = y,
            gam = gam,
            rho = rho,
            mu = mu,
            sig = sig,
            LAM = LAM,
            THETA = THETA,
            OMEGA = OMEGA,
            bet = bet,
            bete = d.btil,
            varphi = d.varphi,
            gr = d.gr,
            sigma_lag_coeff = sigma_lag_coeff
        )
    end

    # Simulate the paths used in the four panels of Figure 2
    function simulate_figure2(ss; T=60, r0=-0.01, rho_r=0.5)
        n = T + 1
        rhat = [r0 * rho_r^t for t in 0:T]

        yhat        = zeros(n)
        muhat       = zeros(n)
        muhat_asset = zeros(n)

        y_next         = 0.0
        mu_next        = 0.0
        mu_asset_next  = 0.0

        # Backward recursion for output and passthrough terms
        for t in n:-1:1
            yhat[t] = ss.THETA * y_next -
                      (1.0 / (ss.gam * ss.y)) * rhat[t] -
                      (ss.LAM / (ss.gam * ss.y)) * mu_next

            muhat[t] = -ss.gam * ss.mu * ss.w * ss.y * (1.0 + ss.gam * ss.rho) * yhat[t] +
                       ss.bete * (mu_next + rhat[t])

            muhat_asset[t] = ss.bete * (mu_asset_next + rhat[t])

            y_next        = yhat[t]
            mu_next       = muhat[t]
            mu_asset_next = muhat_asset[t]
        end

        sigma_total     = zeros(n)
        sigma_asset     = zeros(n)
        sigma_noincome  = zeros(n)

        prev_total    = 0.0
        prev_asset    = 0.0
        prev_noincome = 0.0

        # Build the three versions of inequality dynamics
        for t in 1:n
            sigma_total[t] =
                ss.LAM * muhat[t] - ss.gam * ss.y * (ss.THETA - 1.0) * yhat[t] +
                ss.sigma_lag_coeff * prev_total

            sigma_asset[t] =
                ss.LAM * muhat_asset[t] +
                ss.sigma_lag_coeff * prev_asset

            sigma_noincome[t] =
                ss.LAM * muhat[t] +
                ss.sigma_lag_coeff * prev_noincome

            prev_total    = sigma_total[t]
            prev_asset    = sigma_asset[t]
            prev_noincome = sigma_noincome[t]
        end

        return (
            t = collect(0:T),
            rhat = rhat,
            yhat = yhat,
            muhat = muhat,
            muhat_asset = muhat_asset,
            sigma_total = sigma_total,
            sigma_asset = sigma_asset,
            sigma_noincome = sigma_noincome
        )
    end

    let
        ss_fig2  = steady_state_fig2(baseline_fig2)
        sim_fig2 = simulate_figure2(ss_fig2; T=60, r0=-0.01, rho_r=0.5)

        # The paper shows only the first periods, so we plot the first 11 points
        nplot_fig2 = 11
        tplot_fig2 = sim_fig2.t[1:nplot_fig2]

        fig2_p1 = plot(
            tplot_fig2, 100 .* sim_fig2.rhat[1:nplot_fig2],
            lw = 2.5,
            xlabel = "t",
            ylabel = "% pts",
            title = "a. Real rate path",
            legend = false,
            grid = false
        )

        fig2_p2 = plot(
            tplot_fig2, 100 .* sim_fig2.muhat[1:nplot_fig2],
            lw = 2.5,
            xlabel = "t",
            ylabel = "Ã—100",
            title = "b. Passthrough Î¼Ì‚_t",
            label = "total",
            grid = false
        )
        plot!(
            fig2_p2,
            tplot_fig2, 100 .* sim_fig2.muhat_asset[1:nplot_fig2],
            lw = 2.5,
            ls = :dash,
            label = "asset market only"
        )

        fig2_p3 = plot(
            tplot_fig2, 100 .* sim_fig2.yhat[1:nplot_fig2],
            lw = 2.5,
            xlabel = "t",
            ylabel = "Ã—100",
            title = "c. Output Å·_t",
            legend = false,
            grid = false
        )

        fig2_p4 = plot(
            tplot_fig2, 100 .* sim_fig2.sigma_total[1:nplot_fig2],
            lw = 2.5,
            xlabel = "t",
            ylabel = "Ã—100",
            title = "d. Inequality Î£Ì‚_t",
            label = "full effect",
            grid = false
        )
        plot!(
            fig2_p4,
            tplot_fig2, 100 .* sim_fig2.sigma_asset[1:nplot_fig2],
            lw = 2.5,
            ls = :dash,
            label = "asset market only"
        )
        plot!(
            fig2_p4,
            tplot_fig2, 100 .* sim_fig2.sigma_noincome[1:nplot_fig2],
            lw = 2.5,
            ls = :dot,
            label = "no income-risk channel"
        )

        fig2 = plot(
            fig2_p1, fig2_p2, fig2_p3, fig2_p4,
            layout = (1, 4),
            size = (1500, 340)
        )

        save_figure(fig2, "figure2")

        fig2
    end
end

