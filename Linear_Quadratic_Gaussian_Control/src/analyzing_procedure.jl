#----------------------------------------------------------------------------------#

function analyzing_procedure(
    d::Int64, 
    T::Int64, 
    ρ::Real, 
    N::Vector{Int64};
    )::Tuple{
        <:Real,
        <:Real, 
        Tuple{<:Vector{<:Vector{<:Real}}, <:Vector{<:Vector{<:Real}}, 
            <:Vector{<:Vector{<:Real}}},
        Tuple{<:Vector{<:Vector{<:Real}}, <:Vector{<:Vector{<:Real}}, 
            <:Vector{<:Vector{<:Real}}}
        }

    P = ρ * ones(2 * T + 1)

    # generate linear quadratic gaussian system and store covariance matrices
    system = generate_extended_LQGS(d, T)
    Z_optimal = [system.X, system.W..., system.V...]

    # compute optimal value and optimal robust value
    P_optimal = f(system, Z_optimal)
    _, P_optimal_robust = Frank_Wolfe_Method_LQGS(system, Z_optimal, Z_optimal, P)
    println()

    # initialize output variables
    best_case_trajectory_collection = Vector{Vector{Float64}}(undef, N[2])
    best_case_computation_time_collection = Vector{Vector{Float64}}(undef, N[2])
    best_case_error_collection = Vector{Vector{Float64}}(undef, N[2])

    adaptive_trajectory_collection = Vector{Vector{Float64}}(undef, N[2])
    adaptive_computation_time_collection = Vector{Vector{Float64}}(undef, N[2])
    adaptive_error_collection = Vector{Vector{Float64}}(undef, N[2])

    overall_time = 0

    # start averaging 
    for i = 1:N[2]
        # compute an initial estimate with using N[1] many sample points 
        Z_estimated = initial_estimations(system, N[1])

        Z_iteration_best_case = nothing
        Z_iteration_adaptive = nothing

        # initialize inner output variables
        best_case_trajectory = Vector{Float64}(undef, N[3])
        best_case_computation_time = Vector{Float64}(undef, N[3])
        best_case_error = Vector{Float64}(undef, N[3])

        adaptive_trajectory = Vector{Float64}(undef, N[3])
        adaptive_computation_time = Vector{Float64}(undef, N[3])
        adaptive_error = Vector{Float64}(undef, N[3])


        ns = length(string(N[2]))
        ms = length(string(N[3]))

        # start sampling
        time = @elapsed begin 
            for k = 1:N[3]
                print("\r")
                print("   AVERAGING RUN ", lpad(i, ns, ' '), 
                    " OF $(N[2]) | SAMPLE ITERATION ", lpad(k, ms, ' '),
                    " OF $(N[3]) |", " $(round(100 * (k + (i - 1) * 
                    N[3]) / (N[2] * N[3]), digits = 2))%")

                # assign initial iterate or update estimate
                if k == 1
                    Z_iteration_best_case = Z_estimated
                    Z_iteration_adaptive = Z_estimated
                else
                    Z_estimated = update_estimations(system, Z_estimated, N[1], k)
                end

                # compute inner output variables
                best_case_computation_time[k] = @elapsed Z_iteration_best_case, 
                    best_case_trajectory[k] = Frank_Wolfe_Method_LQGS(system, 
                    Z_iteration_best_case, Z_estimated, P; verbose = false)
                best_case_error[k] = 
                    norm(best_case_trajectory[k] - P_optimal_robust, 2)
            
                adaptive_computation_time[k] = @elapsed Z_iteration_adaptive, _ = 
                    Frank_Wolfe_Method_LQGS_iteration(system, Z_iteration_adaptive, 
                    Z_estimated, P, k)
                adaptive_trajectory[k] = f(system, Z_iteration_adaptive) 
                adaptive_error[k] = 
                    norm(adaptive_trajectory[k] - P_optimal_robust, 2)
            end
        end

        overall_time += time

        println(" (", @sprintf("%9.5f", time), " seconds)")

        # collect output variables
        best_case_trajectory_collection[i] = best_case_trajectory
        best_case_computation_time_collection[i] = best_case_computation_time
        best_case_error_collection[i] = best_case_error

        adaptive_trajectory_collection[i] = adaptive_trajectory
        adaptive_computation_time_collection[i] = adaptive_computation_time
        adaptive_error_collection[i] = adaptive_error
    end

    println("\n   Overall computation time: ", @sprintf("%9.5f", 
        overall_time), " seconds \n   Average computation time: ", 
        @sprintf("%9.5f", overall_time / N_averaging), " seconds")

    best_case_data = (best_case_trajectory_collection, 
        best_case_computation_time_collection, best_case_error_collection)
    adaptive_data = (adaptive_trajectory_collection, 
        adaptive_computation_time_collection, adaptive_error_collection)

    return (P_optimal, P_optimal_robust, best_case_data, adaptive_data)
end

#----------------------------------------------------------------------------------#