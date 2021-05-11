function Compute_deterministic_R0(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   deterministic R0: deterministic squart reproductive number

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    b = par.b[1]
    beta_v = par.beta_v[1]
    gamma = par.gamma[1]

    Rd0 = (beta_p * beta_v * b) / (gamma * (b + r_1) * r_2)

    return Rd0
end
