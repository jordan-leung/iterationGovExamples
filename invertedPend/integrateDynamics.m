function [output] = integrateDynamics(t,x0,controlArgs)

% Get length vectors
N  = length(t);
n  = length(x0);
m = controlArgs.m;

% * Initialize vectors
X      = zeros(n,N);
X(:,1) = x0;

% Get MPC parameters
A = controlArgs.A;
B = controlArgs.B;
Q = controlArgs.Q;
R = controlArgs.R;
P = controlArgs.P;
K = controlArgs.K;
N_MPC = controlArgs.N;
r = controlArgs.r;
xmin = controlArgs.xmin;
xmax= controlArgs.xmax;
umin = controlArgs.umin;
umax = controlArgs.umax;
Gx = controlArgs.Gx;
Gu = controlArgs.Gu;
lambda = controlArgs.lambda;
alpha  = controlArgs.alpha;
options = controlArgs.options;
ws_flag = controlArgs.ws_flag;

% Generate the constant QP matrices (fixed horizon lengths)
[H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);
H_QP = H_MPC;
eigH = eig(H_QP);
options.eigH = eigH;

% Generate minmax vector
zmax = zeros(size(H_MPC,1),1);
zmin = zmax;
for i = 1:N_MPC
    zmin(1 + m*(i-1) : m + m*(i-1)) = umin;
    zmax(1 + m*(i-1) : m + m*(i-1)) = umax;
end
        
% Storage Variables
U = zeros(m,N-1);
numIterHist = zeros(N-1,1);
refHist = zeros(size(r));
VHist = zeros(N-1,1);
U_MPC = zmax*0;

for i=1:N-1
    % Get current state vector
    Xi    = X(:,i);

    % Execute control law at current time-step - If unspecified flag is
    % passed in, the dynamics will integrate with no control law
    ref_i = r(:,i);
    xTilde_i = Xi - Gx*ref_i;
    c_QP = G_MPC*xTilde_i;

   % Run MPC
    if ws_flag
        U_MPC = projGradSolver(H_QP,c_QP,U_MPC,zmin,zmax,options);
    else
        U_MPC = projGradSolver(H_QP,c_QP,U_MPC*0,zmin,zmax,options);
    end


    % Integrate dynamics
    Ui =  U_MPC(1:m);
    X(:,i+1) = A*Xi + B*Ui;
    U(:,i)   = Ui;

    % Store variables
    numIterHist(i) = 1;
    refHist(:,i) = ref_i;
    V_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
    % [V_i,terminalFlag_i] = evaluateCost(Xi,U_MPC,N_MPC,m,A,B,lambda*P,Q,R,Gx*ref_i,Gu*ref_i,alpha,lambda);
    VHist(i) = V_i;
end

% Define outputs
output.X = X;
output.U = U;
output.numIterHist = numIterHist;
output.refHist = refHist;
output.VHist = VHist;
% output.terminalFlag = terminalFlag;
end % end function

function [COST,terminalFlag] = evaluateCost(X0,U,N,m,A,B,P,Q,R,xBar,uBar,alpha,lambda)

% Do the first point outside the loop to initialize
X_i = X0;
U_i = U(1:m);
COST = (X_i - xBar)'*Q*(X_i-xBar) + (U_i - uBar)'*R*(U_i-uBar);
X_i = A*X_i + B*U_i;

% Do 2:N-1
for i = 2:N
    U_i = U(1 + (i-1)*m:m + (i-1)*m);
    COST = COST + (X_i - xBar)'*Q*(X_i-xBar) + (U_i - uBar)'*R*(U_i-uBar);
    X_i = A*X_i + B*U_i;
end

% Do N
terminalCost = (X_i - xBar)'*P*(X_i-xBar);
COST = COST + terminalCost;

% Check terminal condition
if nargin > 11
    if terminalCost < lambda*alpha
        terminalFlag = 1;
    else
        terminalFlag = 0;
    end
else
    terminalFlag = [];
end

end




