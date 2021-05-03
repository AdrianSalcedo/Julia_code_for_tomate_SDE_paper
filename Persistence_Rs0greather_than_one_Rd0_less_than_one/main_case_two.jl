path1 = "/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/"
path2 ="/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/Case_two/"
include(path1 * "Compute_fixed_points.jl")
include(path1 * "Verify_no_extinction_by_noise.jl")
include(path1 * "Case_two/Verify_R0s_greaterthan_R0d_noise_condition.jl")
include(path1 * "Compute_stochastic_R0.jl")
include(path1 * "Case_two/Verify_rd_lessthan_one.jl")
include(path1 * "Verify_rs_persistence_condition.jl")
include(path1 * "Compute_auxiliar_constants_ci.jl")
include(path1 * "Compute_auxiliar_constants_rho_i.jl")
include(path1 * "Verify_persistence_assumption2.jl")



Persistence_parameters = CSV.read(path2 * "Parameter_Persistence_case_two.csv", DataFrame)
Constanst_rho_i = CSV.read(path2 * "Constant_rho_i_case_two.csv", DataFrame)
Constanst_a_i = CSV.read(path2 * "Constants_a_i_case_two.csv", DataFrame)
Constanst_c_i = CSV.read(path2 * "Constants_c_i_case_two.csv", DataFrame)
Endemic_fixed_point = CSV.read(path2 * "endemic_fixed_point_case_two.csv", DataFrame)


R0 = Compute_deterministic_R0(Persistence_parameters)

Rs0 = Compute_stochastic_R0(Persistence_parameters)
Condition_no_extinction_By_noise =
    Verify_no_extinction_by_noise(Persistence_parameters)
    Verify_rd_lessthan_one(Persistence_parameters)
    Verify_R0s_greaterthan_R0d_noise_condition(Persistence_parameters)
    Verify_rs_persistence_condition(Persistence_parameters)
    Verify_persistence_assumption2(
            Constanst_c_i,Constanst_rho_i,Persistence_parameters,
            Endemic_fixed_point)
