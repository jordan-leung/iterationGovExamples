clear
close all
clc

addpath('/Users/jordan/Documents/GitHub/optimizationFunctions');
addpath('/Users/jordan/Documents/GitHub/extraFunctions')
addpath('/Users/jordan/Documents/GitHub/iterationGovernor');

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% Set run flag
RUN_FLAG = 'standard';

% Set warm-start flag
ws_flag = 1;
% ws_flag = 0;

% savedataflag
% saveDataFlag = 0;
saveDataFlag = 1;

% Run common file to set simulation parameters
setSimulationParameters;

% Get ROA data
load('ROAData.mat');
VBar = 0.9*min(costNegVec);


%% Calculate super special constants :) 
[H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);

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
ellStar_cs = (log(gamma-beta)-log(mu*sigma*b))/log(eta)
ellStar_ws = (log(1-beta) - log(mu*rho + omega*(1-beta)))/log(eta)
ellCeil = ceil(ellStar_ws);
z1 = mu/(1 - (eta^ellCeil)*omega);
z2 = (1-beta)/(rho*eta^ellCeil);
zeta = (z1+z2)/2;
epsilon = max([beta + zeta*rho*eta^ellCeil, mu/zeta + omega*eta^ellCeil]);

% Package constants
const.beta = beta;
const.eta = eta;
const.mu = mu;
const.b = b;
const.sigma = sigma;
const.gamma = gamma;
const.omega = omega;
const.rho = rho;
const.ellStar_cs = ellStar_cs;
const.ellStar_ws = ellStar_ws;
const.r_psi = sqrt(VBar);
const.r_phi = (1-beta)*sqrt(VBar)/mu;
const.zeta = zeta;
const.epsilon = epsilon;

% Set optimizer options
if ws_flag
    ellStar = ceil(ellStar_ws);
else
    ellStar = ceil(ellStar_cs);
end

% options.MaxIter = 1000;
options.xTol = 1e-16;
v0 = 1;

% * Add controller arguments to a structure
controlArgs.n = n;
controlArgs.m = m;
controlArgs.A = A;
controlArgs.B = B;
controlArgs.P = P;
controlArgs.Q = Q;
controlArgs.R = R; 
controlArgs.K = K;
controlArgs.N = N_MPC; % Horzion length
controlArgs.T = T_MPC; % Sampling period
controlArgs.r = r;
controlArgs.xmin = xmin;
controlArgs.xmax = xmax;
controlArgs.umin = umin;
controlArgs.umax = umax;
controlArgs.Gx = Gx;
controlArgs.Gu = Gu;
controlArgs.Gr = Gr;
controlArgs.lambda = lambda;
controlArgs.gamma = gamma;
controlArgs.alpha = alpha;
controlArgs.nu_0 = v0;
controlArgs.options = options;
controlArgs.ws_flag = ws_flag;
controlArgs.const = const;

%% Execute the simulation

% Generate grid of initial conditions
NSample = 30;
x1lim = 0.65;
x2lim = 0.2;
x1vec = linspace(-x1lim,x1lim,NSample);
x2vec = linspace(-x2lim,x2lim,NSample);
NTotal = NSample^2;
XMatrix = zeros(n,NTotal);
count = 1;
for i = 1:NSample
    for j = 1:NSample
        XMatrix(:,count) = [x1vec(i); x2vec(j)];
        count = count + 1;
    end
end


% Simulation
batchData.controlArgs = controlArgs;
batchData.t = t;
switch RUN_FLAG
    case 'standard'
        % Run the first leg 
        iOuter = 1;
        X0 = XMatrix(:,iOuter);
        [output_i] = integrateDynamics(t,X0,controlArgs);
        output_i.X0 = X0;

        % Initialize the structure
        batchData.output = repmat(output_i,NTotal,1);

        % Run main iteration loop
        for iOuter = 2:NTotal
            if mod(iOuter,10) == 0
                iOuter
            end
            X0 = XMatrix(:,iOuter);
            [output_i] = integrateDynamics(t,X0,controlArgs);
            output_i.X0 = X0;
            batchData.output(iOuter) = output_i;
        end

        % Save data 
        if saveDataFlag
            savestr = ['batchData_tracking_N',num2str(N_MPC)];
            save(savestr,'batchData');
        end
    case 'ig'
        % Run
        [output] = integrateDynamics_ig(t,X0,controlArgs);
        plotRoutine_ig
        output.controlArgs = controlArgs;
        output.t = t;
        if saveDataFlag
            savestr = ['batchData_ig_N',num2str(N_MPC)];
            save(savestr,'batchData');
        end
end
