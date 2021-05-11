using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using IterableTables, DataFrames, DataTables
using StochasticDiffEq
using Distributions

path1 = "/home/gabrielsalcedo/Github/Julia_code_for_tomate_SDE_paper/"
path2 = "Persistence_Rs0greather_than_one_Rd0_less_than_one/"
path = path1 * path2

Parameters = CSV.read(path * "Parameter_Persistence_case_one.csv", DataFrame)
#=
beta_p = Parameters.beta_p[1]
r_1 =  Parameters.r_1[1]
r_2 =  Parameters.r_2[1]
b =  Parameters.b[1]
beta_v =  Parameters.beta_v[1]
theta =  Parameters.theta[1]
mu =  Parameters.mu[1]
gamma =  Parameters.gamma[1]
N_v =  Parameters.N_v[1]
sigma_L =  Parameters.sigma_L[1]
sigma_I =  Parameters.sigma_I[1]
sigma_v = Parameters.sigma_v[1]
=#

beta_p = Parameters.beta_p[1]
r_1 =  Parameters.r_1[1]
r_2 =  Parameters.r_2[1]
b =  Parameters.b[1]
beta_v =  Parameters.beta_v[1]
theta =  Parameters.theta[1]
mu =  Parameters.mu[1]
gamma =  Parameters.gamma[1]
N_v =  Parameters.N_v[1]
sigma_L =  Parameters.sigma_L[1]
sigma_I =  Parameters.sigma_I[1]
sigma_v = Parameters.sigma_v[1]




u_0 = [100.0, 0.0, 0.0, 3.0, 4.0]
T = 1000.0
time = (0.0, T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range(0.0, T, step = 1.0)


R0 = (beta_p * beta_v * b) / (gamma * (b + r_1) * r_2)
Rs0 = R0 - (1 / 2) * ((sigma_L + sigma_I) ^ 2 -
    sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * gamma))
################################################################################

function F_Det(du,u,p,t)
 @inbounds begin
        du[1] = - beta_p * u[1] * u[5] / N_v + r_1 * u[2] + r_2 * u[3]
        du[2] = beta_p * u[1] * u[5] / N_v - (b + r_1) * u[2]
        du[3] = b * u[2] - r_2 * u[3]
        du[4] = - beta_v * u[4] * u[3] / N_p - gamma * u[4] + (1 - theta) * mu
        du[5] = beta_v * u[4] * u[3] / N_p - gamma * u[5] + theta * mu
    end
    nothing
end

function F_Drift(du,u,p,t)
     @inbounds begin
        du[1] = - beta_p * u[1] * u[5] / N_v + r_1 * u[2] + r_2 * u[3]
        du[2] = beta_p * u[1] * u[5] / N_v - (b + r_1) * u[2]
        du[3] = b * u[2] - r_2 * u[3]
        du[4] = - beta_v * u[4] * u[3] / N_p - gamma * u[4] + (1 - theta) * mu
        du[5] = beta_v * u[4] * u[3] / N_p - gamma * u[5] + theta * mu
    end
    nothing

end

function G_Diffusion(du,u,p,t)
    @inbounds begin
        du[1,1] = u[1] * (sigma_L * u[2] + sigma_I * u[3]) / N_p
        du[1,2] = 0
        du[2,1] = - sigma_L * u[1] * u[2] / N_p
        du[2,2] = 0
        du[3,1] = - sigma_I * u[1] * u[3] / N_p
        du[3,2] = 0
        du[4,1] = 0
        du[4,2] = - sigma_v * u[4]
        du[5,1] = 0
        du[5,2] = - sigma_v * u[5]
    end
    nothing
end

################################################################################
######################### Solution computation #################################
########################## Deterministic SolutionPage 241 of The Threshold Behaviour of Epidemic Models ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det)
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,
    noise_rate_prototype = zeros(5,2))
sol = solve(prob_sde_tomato_sys,EM(),dt = dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "R_s = $Rs0", grid = false, showaxis = false, bottom_margin = - 50Plots.px)
p1 = plot(det_sol,vars = (1),color = "blue")
p1 = plot!(sol,vars = (1),color = "darkgreen",title = "Susc. p.")
p2 = plot(det_sol,vars = (2),color = "blue")
p2 = plot!(sol,vars = (2),color = "darkorange", title = "Lat. p.")
p3 = plot(det_sol,vars = (3),color = "blue")
p3 = plot!(sol,vars = (3),color = "darkred",title = "Infec. p.")
p4 = plot(det_sol,vars = (4),color = "blue")
p4 = plot!(sol,vars = (4),color = "green", title = "Susc. v.")
p5 = plot(det_sol,vars = (5),color = "blue")
p5 = plot!(sol,vars = (5),color = "red",title = "Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([[A B C]; [D E F]]),label = "")

################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################
Datos=DataFrame()
j = 0
trajectories = 1
while j <= 10000
    monte_prob = MonteCarloProblem(prob_sde_tomato_sys)
    sim = solve(monte_prob, SROCKC2(),dt= dt,EnsembleThreads(),trajectories=trajectories)
    component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
    component = transpose(component) #transpose to obtain any*5 data matrix
    component = vcat(component...) #to obtain shape for dataframe
    component = vcat(component...) # again do a reshape
    variables = DataFrame(component) # define first data frame
    Datos_aux = DataFrame(t = t_s, S_p = variables[:,1], I_p = variables[:,3], I_v = variables[:,5]) #only some variables
    Datos = append!(Datos, Datos_aux) #append the data in the loop
    j+=1
    println("acepted =",j)
end
#CSV.write("/home/gabrielsalcedo/Dropbox/Artículos/JuliaPro_code/Noise_Extinction/Trajectories//Data_new.csv",Datos)
#CSV.write("D://Data_Persistence.csv",Datos)
