function Verify_r0_grather_one_condition(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   stochastic R0: verify Rs0 persistence condition

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    r_2 = par.r_2[1]
    beta_v = par.beta_v[1]
    gamma = par.gamma[1]
    theta = par.theta[1]

    R0 = Compute_deterministic_R0(par)
    test = R0 > 1
    return test
end
