function Compute_fixed_points(par)
    # Input:
    #   par: vector with parameters
    # Output:
    #   endemic fixed point: vector with endemic fixed point coordinates

    beta_p = par.beta_p[1]
    r_1 = par.r_1[1]
    b = par.b[1]
    r_2 = par.r_2[1]
    beta_v = par.beta_v[1]
    theta = par.theta[1]
    mu = par.mu[1]
    gamma = par.gamma[1]
    sigma_v = par.sigma_v[1]
    epsilon = par.epsilon[1]
    sigma_L = par.sigma_L[1]
    sigma_I = par.sigma_I[1]
    N_v = par.N_v[1]
    N_p = par.N_p[1]

    a_1 = b * gamma * theta * beta_p + gamma * theta * beta_p * r_2 + b *
     beta_p * beta_v
    a_2 = b * gamma * theta * N_v * r_2 + gamma * theta * N_v * r_1 * r_2
    a_3 = b * gamma * beta_p + b * beta_p * beta_v + gamma * beta_p * r_2
    a_4 = b * gamma * N_v * r_2 + gamma * N_v * r_1 * r_2

    w_1 = (-(a_4-N_v*a_1)+sqrt((a_4-N_v*a_1)^2+4*a_3*N_v*a_2))/(2*a_3)

    S_p_aster = (N_p * N_v * r_2 * (b + r_1)) / (b * beta_p * w_1 + b * N_v *
     r_2 + beta_p * r_2 * w_1 + N_v * r_1 * r_2)
    L_p_aster = (beta_p * N_p * r_2 * w_1) / (b * beta_p * w_1 + b * r_2 *
     N_v + beta_p * r_2 * w_1 + N_v * r_1 * r_2)
    I_p_aster = N_p - (S_p_aster + L_p_aster)
    S_v_aster = ((1 - theta) * mu * (beta_p * b * w_1 + b * r_2 * N_v + beta_p *
        r_2 * w_1 + N_v * r_1 * r_2)) / (b * gamma * beta_p * w_1 + b * gamma *
        N_v * r_2 + b * beta_p * beta_v * w_1 + gamma * beta_p * r_2 * w_1 +
        gamma * N_v * r_1 * r_2)
    I_v_aster = w_1

    endemic_fixed_point = DataFrame(S_p_aster = S_p_aster,
        L_p_aster = L_p_aster, I_p_aster = I_p_aster, S_v_aster = S_v_aster,
        I_v_aster = I_v_aster)
return endemic_fixed_point
end
