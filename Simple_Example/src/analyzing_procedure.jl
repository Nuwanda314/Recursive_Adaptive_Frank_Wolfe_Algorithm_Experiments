#----------------------------------------------------------------------------------#

function start_analysis(
    P::problem_extended,
    N::Vector{Int64}, 
    T::theoretic_data;
    quadratic::Bool = false,
    convex_hull_representation::Bool = false
)::Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}

    best_case_errors = Vector{Vector{Float64}}(undef, N[2])
    adaptive_errors = Vector{Vector{Float64}}(undef, N[2])

    for i = 1:N[2]
        # initialize output variables
        best_case_errors[i] = Vector{Float64}(undef, N[3] + 1)
        adaptive_errors[i] = Vector{Float64}(undef, N[3] + 1)

        A, B = initial_domain_estimation(P, N[1]; 
            convex_hull_representation = convex_hull_representation)

        # in the quadratic convergence setting extend the domain
        if quadratic
            erosion_term = T.c * (2 ./ (2 .+ N[1])) .^ (1 / T.r)

            A += erosion_term
            A -= erosion_term
        end

        x_adaptive = (A + B) / 2 # initial iterate

        best_case_errors[i][1] = abs(optimal_data(P; A, B)[2] - P.P_optimal)
        adaptive_errors[i][1] = abs(f(P, x_adaptive) - P.P_optimal)

        for k = 0:(N[3] - 1)
            if quadratic
                A -= erosion_term
                A += erosion_term
            end

            # depending on the domain approximation style update the new domain 
            # approximation accordingly
            if convex_hull_representation
                A, B = update_domain_estimation(P, A, B; 
                    convex_hull_representation = convex_hull_representation)
            else
                A, B = update_domain_estimation(P, A, B; N = N[1], k = k)
            end

            # in the quadratic convergence setting extend the domain
            if quadratic
                erosion_term = T.c * (2 ./ (2 .+ N[1] + k)) .^ (1 / T.r)

                A += erosion_term
                A -= erosion_term
            end

            best_case_errors[i][k + 2] = 
                abs(optimal_data(P; A, B)[2] - P.P_optimal)

            # adaptive procedure
            s = LMO(P, A, B, x_adaptive)

            λ = 2 / (k + 2)
            x_adaptive = x_adaptive + λ * (s - x_adaptive)

            adaptive_errors[i][k + 2] = abs(f(P, x_adaptive) - P.P_optimal)
        end
    end

    return (best_case_errors, adaptive_errors)
end

#----------------------------------------------------------------------------------#