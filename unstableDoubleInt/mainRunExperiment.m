clear
close all
clc

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

%% * Set simulation parameters

% Set run flag
RUN_FLAG = 'standard'; % run rMPC
% RUN_FLAG = 'ig'; % run IG-MPC
% RUN_FLAG = 'sg'; % run SG-MPC (stability gov. w/ 'optimal' MPC)

% Set warm-start flag
% ws_flag = 1;
ws_flag = 0;

% Set save data flag
% saveDataFlag = 0;
saveDataFlag = 1;

% Run common file to set simulation parameters
setSimulationParameters;

% Set initial condition
icflag = 3; % 1 through 4
switch icflag
    case 1
        X0 = [-0.6; 0];
        eqFlag = 1; % eqflag indicates whether or not the provided IC is an equilibrium of the system 
    case 2
        X0 = [0.7; 0];
        eqFlag = 1;
    case 3 
        X0 = [1; 0.075];
        eqFlag = 0;
    case 4
        X0 = [-0.9; -0.075];
        eqFlag = 0;
    case 5
        X0 = [-1; 0];
        eqFlag = 1;
end

% Get ROA data
load('./Data/ROAData.mat');
VBar = 0.8*min(costNegVec); % get size of ROA as 80% of the smallest value outside Gamma_N


%% Calculate constants
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

%% Execute the simulation

% Set optimizer options
if ws_flag
    ellStar = ceil(ellStar_ws);
else
    ellStar = ceil(ellStar_cs);
end
% options.MaxIter = 1000;
options.xTol = 1e-6;
options.MaxIter = 1e6;
v0 = X0(1);

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
controlArgs.v0 = v0;
controlArgs.options = options;
controlArgs.ws_flag = ws_flag;
controlArgs.const = const;
controlArgs.eqFlag = eqFlag;

% Simulation
switch RUN_FLAG
    case 'standard'
          % Run
        [output] = integrateDynamics(t,X0,controlArgs);
        plotRoutine
        output.controlArgs = controlArgs;
        output.t = t;
        if saveDataFlag
            savestr = ['./Data/output_tracking_N',num2str(N_MPC),'_ic',num2str(icflag)];
            save(savestr,'output');
        end
    case 'ig'
        % Run
        [output] = integrateDynamics_ig(t,X0,controlArgs);
        plotRoutine_ig
        output.controlArgs = controlArgs;
        output.t = t;
        if saveDataFlag
            savestr = ['./Data/output_ig_N',num2str(N_MPC),'_ic',num2str(icflag),'_ws',num2str(ws_flag)];
            save(savestr,'output');
        end
    case 'sg'
        % Run
        [output] = integrateDynamics_sg(t,X0,controlArgs);
        plotRoutine_ig
        output.controlArgs = controlArgs;
        output.t = t;
        if saveDataFlag
            savestr = ['./Data/output_sg_N',num2str(N_MPC),'_ic',num2str(icflag),'_ws',num2str(ws_flag)];
            save(savestr,'output');
        end
end
