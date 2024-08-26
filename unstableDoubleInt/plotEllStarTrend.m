clc
clear all
close all

%% Run main loop

% Run common file to set simulation parameters
setSimulationParameters;
NVec = 1:100;
NSample = length(NVec);

% Initialize
ellStar_ws = zeros(NSample,1);
ellStar_cs = zeros(NSample,1);

for iOuter = 1:NSample
    % Set matrices
    N_MPC = NVec(iOuter);
    [H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,AHat,BHat] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);
    
    % Intermediate matrices
    condH = cond(H_MPC);
    sqrtInvH = sqrtm(inv(H_MPC));
    Xi = [eye(m), zeros(m,m*(N_MPC-1))];
    BBar = B*Xi;

    % Constants
    beta = sqrt(1- weightEig(Q,W_MPC,'-'));
    eta = (condH-1)/(condH+1);
    mu = norm(sqrtm(W_MPC)*BBar);
    b = norm(sqrtInvH);
    if lambda > 1
        sigma = sqrt(weightEig(W_MPC,Q,'+') -1);
    else
        sigma = sqrt(weightEig(W_MPC,P,'+') -1);
    end
    gamma = (1+beta)/2;
    omega = 1 + b*norm(sqrtInvH*G_MPC*BBar);
    rho = b*norm(sqrtInvH*G_MPC*(A-eye(n))*sqrtm(inv(Q))) ...
        + b*sqrt(weightEig(G_MPC*BBar,H_MPC,'+')*sigma^2);
    ellStar_cs(iOuter) = (log(gamma-beta)-log(mu*sigma*b))/log(eta);
    ellStar_ws(iOuter) = (log(1-beta) - log(mu*rho + omega*(1-beta)))/log(eta);
end

semilogy(NVec,ellStar_ws,'r')
hold on; box on; grid on;
semilogy(NVec,ellStar_cs,'b')

% Save the data (will plot in other file)
save('ellStarTrend','NVec','ellStar_ws','ellStar_cs');


