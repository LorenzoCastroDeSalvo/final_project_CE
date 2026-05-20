using Pkg

const PROJECT_ROOT = normpath(@__DIR__)

println("Activating project at: ", PROJECT_ROOT)

Pkg.activate(PROJECT_ROOT)
Pkg.instantiate()

cd(PROJECT_ROOT) do
    mkpath("figures")
    mkpath("output")

    println("Running policy value plots...")
    include(joinpath(PROJECT_ROOT, "scripts", "policy_value_plots.jl"))

    println("Running HANK figures...")
    include(joinpath(PROJECT_ROOT, "scripts", "hank-figures.jl"))
end

println("All scripts completed successfully.")