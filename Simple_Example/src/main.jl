#----------------------------------------------------------------------------------#
# PACKAGES                                                                         #
#----------------------------------------------------------------------------------#

using Pkg
Pkg.add(["CairoMakie", "Random", "Statistics"])

using CairoMakie, Random, Statistics

Random.seed!(42);

#----------------------------------------------------------------------------------#
# INCLUDE FUNCTIONS                                                                #
#----------------------------------------------------------------------------------#

# required packages: -
include("structures.jl");

# required packages: -
include("objective_function.jl");

# required packages: -
include("oracles.jl");

# required packages: Random, Statistics 
include("domain_approximation.jl");

# required packages: - 
include("analyzing_procedure.jl");

# required packages: CairoMakie, Statistics 
include("plot_generation.jl");

#----------------------------------------------------------------------------------#
# PLOT GENERATION                                                                  #
#----------------------------------------------------------------------------------#

a = 0;
b = 1;
x_star = 2;

P = extend_problem(problem(a, b, x_star));

m = 20;
τ = 2;

T = generate_theoretic_data(P, m, τ; convex_hull_representation = false);

N_domain = 10;
N_averaging = 25;
N_sampling = 10000;

N = [N_domain, N_averaging, N_sampling];

(best_case_errors, adaptive_errors) = start_analysis(P, N, T; quadratic = false, 
    convex_hull_representation = false);

#----------------------------------------------------------------------------------#

pa = generate_error_comparison_plot(
    N, 
    T,
    [best_case_errors, adaptive_errors];
    y_lim_lower = 10^-4, 
    y_lim_upper = 10^4,
    title = "(a) Error Comparison - Moment-Based Domain Approximation"
)

save(joinpath(@__DIR__, "..", "figures", "EX1_(a).png"), pa);

#----------------------------------------------------------------------------------#

T = generate_theoretic_data(P, m, τ; convex_hull_representation = true);

(best_case_errors, adaptive_errors) = start_analysis(P, N, T; quadratic = false, 
    convex_hull_representation = true);

#----------------------------------------------------------------------------------#

pb = generate_error_comparison_plot(
    N, 
    T,
    [best_case_errors, adaptive_errors];
    y_lim_lower = 10^-4, 
    y_lim_upper = 10^4,
    title = "(b) Error Comparison - Convex Hull Domain Approximation"
)

save(joinpath(@__DIR__, "..", "figures", "EX1_(b).png"), pa);

#----------------------------------------------------------------------------------#

x_star = 0.5;

P = extend_problem(problem(a, b, x_star));
T = generate_theoretic_data(P, m, τ; convex_hull_representation = true);

N_sampling = 1000000;

N = [N_domain, N_averaging, N_sampling];

(best_case_errors, adaptive_errors) = start_analysis(P, N, T; quadratic = true, 
    convex_hull_representation = true);

#----------------------------------------------------------------------------------#

pc = generate_accelerated_error_plot(
    N, 
    T,
    adaptive_errors;
    spacing = 10,
    y_lim_lower = 10^-14, 
    y_lim_upper = 10^0,
    title = "(c) Accelerated Convergence with Optimal Solution in the Interior"
)

save(joinpath(@__DIR__, "..", "figures", "EX1_(c).png"), pa);

#----------------------------------------------------------------------------------#

pd = generate_accelerated_error_plot(
    N, 
    T,
    nothing;
    spacing = 10,
    y_lim_lower = 10^-6, 
    y_lim_upper = 10^0,
    title = "(d) Visualization of the Improved Theoretical Bound"
)

save(joinpath(@__DIR__, "..", "figures", "EX1_(d).png"), pa);

#----------------------------------------------------------------------------------#