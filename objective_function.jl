#----------------------------------------------------------------------------------#

function f(
    system::LQGS_extended, 
    Z::Vector{<:Matrix{<:Real}} # Z = [X, W[1],..., W[T], V[1],..., V[T]]
    )::Real

    J = 0
    Σ_predicted = Z[1] # initialize covariance prediction

    # TODO: DESCRIPTION OF LOOP
    for t = 1:system.T
        c₁ = system.C[t] * Σ_predicted

        # compute inverse via \ operator
        c₂ = inv(c₁ * system.C[t]' + Z[t + system.T + 1])
        Σ = Σ_predicted - c₁' * c₂ * c₁

        J += tr((system.Q[t] - system.P[t]) * Σ) + tr(system.P[t] * Σ_predicted)

        Σ_predicted = system.A[t] * Σ * system.A[t]' + Z[t + 1]
    end

    # add last term, where Q[T + 1] - P[T + 1] = 0 by construction
    J += tr(system.P[system.T + 1] * Σ_predicted) # default return
end

#----------------------------------------------------------------------------------#

function ∇f(
    system::LQGS_extended, 
    Z::Vector{<:Matrix{<:Real}}, # Z = [X, W[1], ..., W[T], V[1], ..., V[T]]
    k::Int64                     # k ∈ [1, ..., 2T + 1]
    )::Matrix{<:Real}

    # depending on the input index k, distinguish between the differentiation with 
    # respect to the variable X, W[1], ..., W[T], V[1], ..., V[T] which are the 
    # components of the input variable Z. Use package 'ForwardDiff' to differentite.

    # Z[k] = X
    if k == 1
        function ∇f_X(
            X::Matrix{<:Real}      # X = Z[1]
            )::Real

            J = 0
            Σ_predicted = X # initialize covariance prediction

            for t = 1:system.T
                c₁ = system.C[t] * Σ_predicted

                # compute inverse via \ operator
                c₂ = (c₁ * system.C[t]' + Z[t + system.T + 1]) \ I
                Σ = Σ_predicted - c₁' * c₂ * c₁

                J += tr((system.Q[t] - system.P[t]) * Σ) + 
                    tr(system.P[t] * Σ_predicted)

                Σ_predicted = system.A[t] * Σ * system.A[t]' + Z[t + 1]
            end

            # add last term, where Q[T + 1] - P[T + 1] = 0 by construction
            J += tr(system.P[system.T + 1] * Σ_predicted)
        end

        return ForwardDiff.gradient(∇f_X, Z[k])

    # Z[k] = W[1], ..., W[T]
    elseif k ∈ 2:(system.T + 1)
        function ∇f_W(
            W::Matrix{<:Real}      # W = Z[k]
            )::Real

            J = 0 
            Σ_predicted = Z[1] # initialize covariance prediction

            for t = 1:system.T
                c₁ = system.C[t] * Σ_predicted

                # compute inverse via \ operator
                c₂ = (c₁ * system.C[t]' + Z[t + system.T + 1]) \ I
                Σ = Σ_predicted - c₁' * c₂ * c₁

                J += tr((system.Q[t] - system.P[t]) * Σ) + 
                    tr(system.P[t] * Σ_predicted)

                # if k = t + 1, replace Z[k] with W
                Σ_predicted = system.A[t] * Σ * system.A[t]' + 
                    (k == t + 1 ? W : Z[t + 1])
            end

            # add last term, where Q[T + 1] - P[T + 1] = 0 by construction
            J += tr(system.P[system.T + 1] * Σ_predicted)
        end

        return ForwardDiff.gradient(∇f_W, Z[k])

    # Z[k] = V[1], ..., V[T]
    elseif k ∈ (system.T + 2):(2 * system.T + 2)
        function ∇f_V(
            V::Matrix{<:Real}      # V = Z[k]
            )::Real

            J = 0
            Σ_predicted = Z[1] # initialize covariance prediction

            for t = 1:system.T
                c₁ = system.C[t] * Σ_predicted

                # compute inverse via \ operator 
                # if k = t + T + 1, replace Z[k] with V 
                c₂ = (c₁ * system.C[t]' + 
                    (k == t + system.T + 1 ? V : Z[t + system.T + 1])) \ I
                Σ = Σ_predicted - c₁' * c₂ * c₁

                J += tr((system.Q[t] - system.P[t]) * Σ) + 
                    tr(system.P[t] * Σ_predicted)

                Σ_predicted = system.A[t] * Σ * system.A[t]' + Z[t + 1]
            end

            # add last term, where Q[T + 1] - P[T + 1] = 0 by construction
            J += tr(system.P[system.T + 1] * Σ_predicted)
        end

        return ForwardDiff.gradient(∇f_V, Z[k])
    end
end

#----------------------------------------------------------------------------------#