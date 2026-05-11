#----------------------------------------------------------------------------------#
# PACKAGES                                                                         #
#----------------------------------------------------------------------------------#

using Pkg   
Pkg.add(["CairoMakie", "Distributions", "ForwardDiff", "JuMP", "LinearAlgebra", 
    "MosekTools", "Printf", "Random", "Statistics"])

using CairoMakie, Distributions, ForwardDiff, JuMP, LinearAlgebra, MosekTools, 
    Printf, Random, Statistics

Random.seed!(42);

#----------------------------------------------------------------------------------#
# INCLUDE OTHER FUNCTIONS                                                          #
#----------------------------------------------------------------------------------#

# required packages: LinearAlgebra
include("structures.jl");

# required packages: LinearAlgebra, Random
include("generate_LQGS.jl");

# required packages: ForwardDiff, LinearAlgebra
include("objective_function.jl");

# required packages: JuMP, LinearAlgebra, MosekTools, Printf
include("Frank_Wolfe_method_LQGS.jl");

# required packages: Distributions, LinearAlgebra, Random
include("gaussian_noise_estimation.jl");

# required packages: LinearAlgebra, Printf, Random
include("analyzing_procedure.jl");

# required packages: CairoMakie, Statistics
include("plot_generation.jl");

#----------------------------------------------------------------------------------#
# PLOT GENERATION                                                                  #
#----------------------------------------------------------------------------------#

d = 10;
T = 10;
ρ = 0.1;

N_initial = 10;
N_averaging = 10;
N_sampling = 1000;

N = [N_initial, N_averaging, N_sampling];

P_optimal, P_optimal_robist, best_case_data, adaptive_data = 
    analyzing_procedure(d, T, ρ, N);

#----------------------------------------------------------------------------------#

pa = generate_trajectory_comparison_plot(
    [P_optimal, P_optimal_robist], 
    [best_case_data[1], adaptive_data[1]], 
    N;
    y_lim_lower = 2.5,
    y_lim_upper = 3.5
)

save(joinpath(@__DIR__, "..", "figures", "EX2_(a).png"), pa);

#----------------------------------------------------------------------------------#

pb = generate_error_comparison_plot(
    [best_case_data[3], adaptive_data[3]], 
    N;
    y_lim_lower = 10^-1.5,
    y_lim_upper = 10^0
)

save(joinpath(@__DIR__, "..", "figures", "EX2_(b).png"), pa);

#----------------------------------------------------------------------------------#
# RUNTIME COMPARISON                                                               #
#----------------------------------------------------------------------------------#

ρ = 0.1;

N_initial = 10;
N_averaging = 10;
N_sampling = 100;

N = [N_initial, N_averaging, N_sampling];

for d ∈ [5, 10, 15]
    for T ∈ [5, 10, 15]
        println("\n   RUNTIME COMPARISON : d = $d, T = $T\n")

        _, _, best_case_data, adaptive_data = analyzing_procedure(d, T, ρ, N)

        best_case_mean = map(x -> mean(x), eachcol(hcat(best_case_data[2]...)'))
        adaptive_mean = map(x -> mean(x), eachcol(hcat(adaptive_data[2]...)'))

        print("\n      MEAN ± SD :", 
            @sprintf("%8.5f", mean(best_case_mean)), " [s] ±")
        print(@sprintf("%8.5f", std(best_case_mean)), " [s] |")
        print(" MIN :", @sprintf("%8.5f", minimum(best_case_mean)), " [s]")

        println()
        print("      MEAN ± SD :", @sprintf("%8.5f", mean(adaptive_mean)), " [s] ±")
        print(@sprintf("%8.5f", std(adaptive_mean)), " [s] |")
        println(" MIN :", @sprintf("%8.5f", minimum(adaptive_mean)), " [s]\n")
    end
end

#----------------------------------------------------------------------------------#