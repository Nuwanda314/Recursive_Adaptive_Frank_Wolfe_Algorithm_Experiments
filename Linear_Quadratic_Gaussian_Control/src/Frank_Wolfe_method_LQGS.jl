#----------------------------------------------------------------------------------#

function SDP_oracle(
    C::Matrix{<:Real}, 
    Z_center::Matrix{<:Real}, 
    ρ::Real
    )::Matrix{<:Real}

    n = size(C)[1]

    model = Model(Mosek.Optimizer)
    set_silent(model)

    # initialize problem variable X and auxilary variable S
    @variable(model, X[1:n, 1:n], PSD)
    @variable(model, S[1:n, 1:n], PSD)

    λ = minimum(eigvals(Z_center))
    @constraint(model, X - λ * I(n) in PSDCone())

    # construct auxilary variable with block structure
    @variable(model, M[1:(2 * n), 1:(2 * n)], Symmetric)
    @constraint(model, M[1:n, 1:n] .== sqrt(Z_center) * X * sqrt(Z_center))
    @constraint(model, M[1:n, (n + 1):end] .== S)
    @constraint(model, M[(n + 1):end, 1:n] .== S)
    @constraint(model, M[(n + 1):end, (n + 1):end] .== I(n))
    @constraint(model, M in PSDCone())

    @constraint(model, tr(X) - 2 * tr(S) ≤ ρ ^ 2 - tr(Z_center))

    # set linear objective for cost matrix C
    @objective(model, Max, tr(C' * X))

    optimize!(model)

    return value.(X)
end

#----------------------------------------------------------------------------------#

function Frank_Wolfe_Method_LQGS(
    system::LQGS_extended, 
    Z_initial::Vector{<:Matrix{<:Real}}, 
    Z_center::Vector{<:Matrix{<:Real}}, 
    ρ::Vector{<:Real}; 
    τ::Real = 1e-3, 
    maximal_iterations::Int64 = 1000, 
    verbose::Bool = true
    )::Tuple{Vector{<:Matrix{<:Real}}, <:Real}

    Z_iteration = copy(Z_initial)
    optimal_value = f(system, Z_iteration)

    if verbose
        print("+----------------------------------------------------------------")
        println("-----------------+")
        print("|   FRANK-WOLFE ALGORITHM STARTED                                ")
        println("                 |")
        print("+----------------------------------------------------------------")
        println("-----------------+")
        print("|   ITERATION ", lpad(0, 4, ' '),
            "  |  OPTIMAL VALUE = ", @sprintf("%14.8f", optimal_value),
            "  |")
    end

    time = @elapsed begin
        for k = 1:(maximal_iterations - 1)
            Z_iteration, frank_wolfe_gap = Frank_Wolfe_Method_LQGS_iteration(system,
                Z_iteration, Z_center, ρ, k)

            optimal_value = f(system, Z_iteration)

            if verbose
                println("  FW GAP = ", @sprintf("%12.8f", frank_wolfe_gap), "   |")
            end

            if frank_wolfe_gap < τ
                break
            end

            if verbose
                print("|   ITERATION ", lpad(k, 4, ' '),
                    "  |  OPTIMAL VALUE = ", @sprintf("%14.8f", optimal_value),
                    "  |")
            end
        end    
    end

    if verbose
        s = @sprintf("%.8f", time)
        n = length(s)
        print("+----------------------------------------------------------------")
        println("-----------------+")
        println("|   FRANK-WOLFE ALGORITHM TERMINATED (", s, " SECONDS)",
            lpad("", 34 - n, ' '), " |")
        print("+----------------------------------------------------------------")
        println("-----------------+")
    end

    return (Z_iteration, optimal_value)
end

#--------------------------------------------------------------------------------#

function Frank_Wolfe_Method_LQGS_iteration(
    system::LQGS_extended, 
    Z_iteration::Vector{Matrix{Float64}}, 
    Z_center::Vector{Matrix{Float64}}, 
    ρ::Vector{<:Real}, 
    k::Int64
    )::Tuple{Vector{<:Matrix{<:Real}}, <:Real}

    frank_wolfe_gap = 0

    X = Vector{Matrix{Float64}}(undef, 2 * system.T + 1)
    for i = 1:(2 * system.T + 1)
        C = ∇f(system, Z_iteration, i)
        X[i] = SDP_oracle(C, Z_center[i], ρ[i])

        frank_wolfe_gap += tr(C' * (X[i] - Z_iteration[i]))
    end
    
    return (Z_iteration + 2 * (X - Z_iteration) / (k + 1), frank_wolfe_gap)
end

#----------------------------------------------------------------------------------#