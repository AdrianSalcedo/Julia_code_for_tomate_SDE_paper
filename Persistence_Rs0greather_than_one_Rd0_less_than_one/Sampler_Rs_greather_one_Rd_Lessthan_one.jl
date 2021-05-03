using CSV
using IterableTables, DataFrames, DataTables
using Distributions
path = "/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/"
include(path * "Compute_fixed_points.jl")
include(path * "Verify_no_extinction_by_noise.jl")
include(path * "Case_two/Verify_R0s_greaterthan_R0d_noise_condition.jl")
include(path * "Compute_stochastic_R0.jl")
include(path * "Case_two/Verify_rd_lessthan_one.jl")
include(path * "Verify_rs_persistence_condition.jl")
include(path * "Compute_auxiliar_constants_ci.jl")
include(path * "Compute_auxiliar_constants_rho_i.jl")
include(path * "Verify_persistence_assumption2.jl")

function Sampler_persistence_parameters(N_p)
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
        N_v = mu/gamma

        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
          sigma_L = sigma_L, sigma_I = sigma_I, sigma_v = sigma_v, N_v = N_v,
           N_p = N_p, epsilon = epsilon)

        par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
         beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
          sigma_L = sigma_L, sigma_I = sigma_I, sigma_v = sigma_v, N_v = N_v,
           N_p = N_p, epsilon = epsilon)

        endemic_fixed_point = Compute_fixed_points(par)
        cond0 = Verify_no_extinction_by_noise(par)
        cond1 = Verify_rd_lessthan_one(par)
        cond2 = Verify_R0s_greaterthan_R0d_noise_condition(par)
        cond3 = Verify_rs_persistence_condition(par)
        auxiliar_constants_c_i =
            Compute_auxiliar_constants_ci(par,endemic_fixed_point)
        auxiliar_constants_rho_i = Compute_auxiliar_constants_rho_i(par)
        cond4, auxiliar_constants_a_i = Verify_persistence_assumption2(
                auxiliar_constants_c_i,auxiliar_constants_rho_i,par,
                endemic_fixed_point
            )
        Test1 = cond0 && cond1 && cond2 && cond3 && cond4
        if Test1 == true
            return par, auxiliar_constants_rho_i, auxiliar_constants_a_i,
            auxiliar_constants_c_i, endemic_fixed_point
        end
    end
end



par, auxiliar_constants_rho_i, auxiliar_constants_a_i,
    auxiliar_constants_c_i, endemic_fixed_point =
        Sampler_persistence_parameters(100)

CSV.write(path * "Case_two//Parameter_Persistence_case_two.csv", par)
CSV.write(path * "Case_two//Constant_rho_i_case_two.csv", auxiliar_constants_rho_i)
CSV.write(path * "Case_two//Constants_a_i_case_two.csv", auxiliar_constants_a_i)
CSV.write(path * "Case_two//Constants_c_i_case_two.csv",auxiliar_constants_c_i)
CSV.write(path * "Case_two//endemic_fixed_point_case_two.csv",endemic_fixed_point)