using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using DataFrames
using StochasticDiffEq

path1 = "/home/gabrielsalcedo/Github/"
path2 = "Julia_code_for_tomate_SDE_paper/Extinction_Noise/"
path = path1 * path2

#include(path * "Dynamics.jl")

beta_p = 0.06
r_1 = 0.92
r_2 = 0.005
b = 0.9
beta_v = 0.2
theta = 0.01
mu = 0.3
gamma = 0.7
sigma_L = 0.7
sigma_I = 0.7
sigma_v = 0.86
N_v = mu/gamma

u_0 = [1.0,0.0,0.0,3.0,3.0]
T = 100.0
time = (0.0,T)
N_p = u_0[1]+u_0[2]+u_0[3]
dt=0.01

tol = 1e-06

################################################################################

cond_1 = beta_p^2/(2*sigma_L^2)+ r_2^2/(2*sigma_I^2)+2*beta_p-r_1
cond_2 = beta_v^2/(2*sigma_v^2)+beta_v-gamma+theta*mu

R0 = sqrt((beta_p*beta_v*b)/(gamma*(b+r_1)*r_2))

println("cond_1=",cond_1);

println("cond_2=",cond_2);

println("R0=",R0);

function F_Det(du,u,p,t)
 @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing
end

function F_Drift(du,u,p,t)
     @inbounds begin
        du[1] = -beta_p*u[1]*u[5]/N_v+r_1*u[2]+r_2*u[3]
        du[2] = beta_p*u[1]*u[5]/N_v-b*u[2]-r_1*u[2]
        du[3] = b*u[2]-r_2*u[3]
        du[4] = -beta_v*u[4]*u[3]/N_p-gamma*u[4]+(1-theta)*mu
        du[5] = beta_v*u[4]*u[3]/N_p-gamma*u[5]+theta*mu
    end
    nothing

end

function G_Diffusion(du,u,p,t)
    @inbounds begin
        du[1,1] = (sigma_L*u[2]*u[1]+sigma_I*u[1]*u[3])/N_p
        du[1,2] = 0
        du[2,1] = -sigma_L*u[1]*u[2]/N_p
        du[2,2] = 0
        du[3,1] = -sigma_I*u[1]*u[3]/N_p
        du[3,2] = 0
        du[4,1] = 0
        du[4,2] = -sigma_v*u[4]
        du[5,1] = 0
        du[5,2] = -sigma_v*u[5]
    end
    nothing
end

######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Det,u_0,time)
det_sol = solve(prob_det)
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift,G_Diffusion,u_0,time,
noise_rate_prototype=zeros(5,2))
stc_sol = solve(prob_sde_tomato_sys,SROCKC2(),dt= dt)
################################################################################
############################ PLot variables ####################################

title = plot(title = "Noise Condition", grid = false, showaxis = false, bottom_margin = -50Plots.px)
p1=plot(det_sol,vars=(1),color="blue")
p1=plot!(stc_sol,vars=(1),color="darkgreen",title="Susc. p.")
p2=plot(det_sol,vars=(2),color="blue")
p2=plot!(stc_sol,vars=(2),color="darkorange", title ="Lat. p.")
p3=plot(det_sol,vars=(3),color="blue")
p3=plot!(stc_sol,vars=(3),color="darkred",title = "Infec. p.")
p4=plot(det_sol,vars=(4),color="blue")
p4=plot!(stc_sol,vars=(4),color="green", title = "Susc. v.")
p5=plot(det_sol,vars=(5),color="blue")
p5=plot!(stc_sol,vars=(5),color="red",title ="Infec. v.")

plot(p1,p2,p3,p4,p5,title,layout = @layout([[A B C]; [D E F]]),label = "")


################################################################################
######################    data  Det Solution    ################################
################################################################################

det_Time = det_sol.t
det_xu = det_sol.u
det_xu_glued = hcat(det_xu...)
Xu2 = det_xu_glued[2:5:end]/N_p
Xu1 = det_xu_glued[1:5:end]/N_p
Xu3 = det_xu_glued[3:5:end]/N_p
Xu4 = det_xu_glued[4:5:end]
Xu5 = det_xu_glued[5:5:end]

det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4, I_v = Xu5)
det_DF1_red = det_DF1#[1:10:end,1:end]
CSV.write(path * "Det_solution_noise_extinction.csv",det_DF1_red)

################################################################################
######################    data  Sto Solution    ################################
################################################################################
stc_Time = stc_sol.t
det_yu = stc_sol.u
det_yu_glued = hcat(det_yu...)
Yu1 = det_yu_glued[1:5:end]/N_p
Yu2 = det_yu_glued[2:5:end]/N_p
Yu3 = det_yu_glued[3:5:end]/N_p
Yu4 = det_yu_glued[4:5:end]
Yu5 = det_yu_glued[5:5:end]

stc_DF1 = DataFrame(t = stc_Time, S_p = Yu1, L_p = Yu2, I_p = Yu3, S_v = Yu4, I_v = Yu5)
stc_DF1_red = stc_DF1[1:25:end,1:end]
CSV.write(path * "Stc_solution_noise_extinction.csv",stc_DF1_red)

################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################

Datos=DataFrame()
j = 0
i = 1
trajectories = 1
t_s = range(0.0,T, step=1.0)
while j <= 10000
    #ensembleprob = EnsembleProblem(prob_sde_tomato_sys)
    #sim = solve(ensembleprob, SROCK1(), dt = dt, abstol=1e-06,dtmin =0.0001,trajectories=trajectories)
    monte_prob = MonteCarloProblem(prob_sde_tomato_sys)
    sim = solve(monte_prob, SROCKC2(),dt= dt,EnsembleThreads(),trajectories=trajectories)
    #plot(sim,idxs=(3))
    component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
    component = transpose(component) #transpose to obtain any*5 data matrix
    component = vcat(component...) #to obtain shape for dataframe
    component = vcat(component...) # again do a reshape
    variables = DataFrame(component) # define first data frame
    Datos_aux = DataFrame(t = t_s, S_p = variables[:,1]/N_p,
        I_p = variables[:,3]/N_p, I_v = variables[:,5]) #only some variables
    #condition1 =  variables[:,1]+variables[:,2]+variables[:,3]
    #condition2 = condition1-N_p*ones(size(condition1))
    #condition3 = maximum(abs.(condition2))
    #condition4 = condition3<=tol
    condition5 = maximum(abs.(variables[end-25:end,3]))
    condition6 = condition5<0.07
    if (condition6 == true)
        Datos = append!(Datos, Datos_aux) #append the data in the loop
        j+=1
        println("acepted =",j)
    end
    println(i)
    i+=1
end

CSV.write(path * "Data_noise_may_10.csv",Datos)
