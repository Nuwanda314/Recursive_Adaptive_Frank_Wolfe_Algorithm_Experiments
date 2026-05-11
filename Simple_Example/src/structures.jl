#----------------------------------------------------------------------------------#

struct problem
    a::Real
    b::Real

    x_star::Real
end

#----------------------------------------------------------------------------------#

struct problem_extended
    a::Real
    b::Real
    
    x_star::Real

    x_optimal::Real
    P_optimal::Real
end

#----------------------------------------------------------------------------------#

function extend_problem(
    P::problem
)::problem_extended

    (x_optimal, P_optimal) = optimal_data(P)

    return problem_extended(P.a, P.b, P.x_star, x_optimal, P_optimal)
end

#----------------------------------------------------------------------------------#

struct theoretic_data
    r::Real

    c::Real
    A::Real
    A_plus::Real
    B::Real

    N_quadratic::Real
end

#----------------------------------------------------------------------------------#

function generate_theoretic_data(
    P::Union{problem, problem_extended},
    m::Int64, 
    τ::Real;
    convex_hull_representation::Bool = false
)::theoretic_data

    if convex_hull_representation
        r = (m - 1) / m
        c = Float64((P.b - P.a) * (2 * 10 ^ τ * factorial(big(m))) ^ (1 / m))
    else
        r = (m - 1) / (2 * m)
        c = Float64(6 * (P.b - P.a) * 
            (2 + (8 * 10 ^ τ * factorial(big(m))) ^ (1 / (2 * m))) ^ 2)
    end

    LL = 6 * (abs(P.b) + abs(P.a)) + 2 * abs(P.x_star) 
    CC = 12 * (abs(P.b) + abs(P.a))

    L_f = 2 * max(abs(P.a - P.x_star), abs(P.b - P.x_star))
    C_f = 2 * (P.b - P.a)

    if convex_hull_representation
        A = 2 * c * L_f + C_f
        A_plus = A

        c_extended = 3 * c * L_f + C_f
    else
        A = 2 * c * LL + CC
        A_plus = 4 * c * L_f + C_f

        c_extended = (3 / 2) * c * L_f + (C_f / 2)
    end

    ρ = min(abs(P.a - P.x_star), abs(P.b - P.x_star))

    B = ((4 * (b - a) / ρ) * 2 * c_extended)^2 + (7 / 2) * c * L_f + (3/ 2) * C_f

    if convex_hull_representation
        N_proposition = ceil(Int64, 2 * (B / c_extended) ^ (1 / r) - 2)
  
    else
        N_proposition = 0
    end
    
    N_quadratic = max(ceil(Int64, 2 * ((6 * c) / ρ) ^ (1 / r) - 2), N_proposition)

    return theoretic_data(r, c, A, A_plus, B, N_quadratic)
end

#----------------------------------------------------------------------------------#