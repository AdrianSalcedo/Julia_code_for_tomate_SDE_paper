function F_Det(du,u,p,t)
 @inbounds begin
        du[1] = - beta_p * u[1] * u[5] + r_1 * u[2] + r_2 * u[3]
        du[2] = beta_p * u[1] * u[5] - b * u[2] - r_1 * u[2]
        du[3] = b * u[2] - r_2 * u[3]
        du[4] = -beta_v * u[4] * u[3] - gamma * u[4] + (1 - theta) * mu
        du[5] = beta_v * u[4] * u[3] - gamma * u[5] + theta * mu
    end
    nothing
end

function F_Drift(du,u,p,t)
     @inbounds begin
         du[1] = - beta_p * u[1] * u[5] + r_1 * u[2] + r_2 * u[3]
         du[2] = beta_p * u[1] * u[5] - b * u[2] - r_1 * u[2]
         du[3] = b * u[2] - r_2 * u[3]
         du[4] = - beta_v * u[4] * u[3] - gamma * u[4] + (1 - theta) * mu
         du[5] = beta_v * u[4] * u[3] - gamma * u[5] + theta * mu
    end
    nothing
end

function G_Diffusion(du,u,p,t)
    @inbounds begin
        du[1,1] = sigma_L * u[2] * u[1] + sigma_I * u[1] * u[3]
        du[1,2] = 0
        du[2,1] = - sigma_L * u[1] * u[2]
        du[2,2] = 0
        du[3,1] = - sigma_I * u[1] * u[3]
        du[3,2] = 0
        du[4,1] = 0
        du[4,2] = - sigma_v * u[4]
        du[5,1] = 0
        du[5,2] = - sigma_v * u[5]
    end
    nothing
end
