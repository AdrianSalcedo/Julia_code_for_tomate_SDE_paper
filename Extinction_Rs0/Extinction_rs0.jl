using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using IterableTables, DataFrames, DataTables
using StochasticDiffEq
using Distributions

path =
    "/home/gabrielsalcedo/Github/Julia_code_for_tomate_SDE_paper/Extinction_Rs0/"

include(path * "Dynamics.jl")

Parameters = CSV.read(path * "Parameter_Extinction_rs0.csv", DataFrame)

beta_p = Parameters.beta_p[1]-0.01
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

u_0 = [ 1.0, 0.0, 0.0, 0.03, 0.01]
    T = 1000.0
time = (0.0, T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range( 0.0, T, step = 1.0)

R0 = (beta_p * beta_v * b) / (gamma * (b + r_1) * r_2)
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * gamma)
        )

################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det,Tsit5(),dt = dt)
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,
noise_rate_prototype=zeros(5,2))
sol = solve(prob_sde_tomato_sys,SROCKC2(),dt= dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(det_sol,vars=(1),color="blue")
p1=plot!(sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(det_sol,vars=(2),color="blue")
p2=plot!(sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(det_sol,vars=(3),color="blue")
p3=plot!(sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(det_sol,vars=(4),color="blue")
p4=plot!(sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(det_sol,vars=(5),color="blue")
p5=plot!(sol,vars=(5),color="red",title ="Infec. v.")

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
