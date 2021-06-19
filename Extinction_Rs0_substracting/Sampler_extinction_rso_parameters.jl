using CSV
using IterableTables, DataFrames, DataTables
using Distributions
path1 =
    "/home/gabrielsalcedo/Github/Julia_code_for_tomate_SDE_paper"
path2 = "/Extinction_Rs0_substracting/"
path = path1 * path2

include(path * "Verify_extinction_by_rso.jl")
include(path * "Verify_R0s_greatherthan_R0d_noise_condition.jl")
include(path * "Compute_stochastic_R0.jl")
include(path * "Verify_rs_extinction_condition.jl")
include(path * "Compute_deterministic_R0.jl")
include(path * "Verify_r0_grather_one_condition.jl")


function Sampler_extinction_rso_parameters(N_p)
    # Input:
    #   N_p: plant size
    # Output:
    #   parameters: vector with parameters that satisfies persistence conditions

    Test1 = false
    while Test1 == false
        beta_p = rand(Uniform(0,1))
        r_1 = rand(Uniform(0,1))
        b = rand(Uniform(0,1))
        r_2 = rand(Uniform((1/2)*(-1+sqrt(2))*b,1))
        beta_v = rand(Uniform(0,1))
        theta = rand(Uniform(0,1))
        mu = rand(Uniform(0,1))
        gamma = rand(Uniform(0,1))
        sigma_v = rand(Uniform(0,1))
        epsilon = rand(Uniform(0,1))
        sigma_L = rand(Uniform(0,1))
        sigma_I = epsilon*sigma_L
        N_p = N_p
        N_v = mu / gamma
        #N_v = N_v / N_v
        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
          sigma_L = sigma_L, sigma_I = sigma_I, sigma_v = sigma_v, N_v = N_v,
           N_p = N_p, epsilon = epsilon)

        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
          sigma_L = sigma_L, sigma_I = sigma_I, sigma_v = sigma_v, N_v = N_v,
           N_p = N_p, epsilon = epsilon)

        cond0 = Verify_extinction_by_rso(par)
        cond1 = Verify_R0s_greatherthan_R0d_noise_condition(par)
        cond2 = Verify_rs_extinction_condition(par)
        cond3 = Verify_r0_grather_one_condition(par)
        Test1 = cond0 && cond1 && cond2 && cond3
        if Test1 == true
            return par
        end
    end
end

par =  Sampler_extinction_rso_parameters(1)

CSV.write(path * "Parameter_Extinction_rs0_substracting.csv", par)
