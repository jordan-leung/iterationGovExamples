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
const = controlArgs.const;
v0 = controlArgs.v0;

% Generate the constant QP matrices (fixed horizon lengths)
[H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,AHat,BHat] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);
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
refHist = zeros(size(r));
kappaHist = zeros(N-1,1);
ellHist = zeros(N-1,1);
ellBarHist = zeros(N-1,1);
errorHist = zeros(N-1,1);
costHist = zeros(N-1,1);
VHist = zeros(N-1,1);
U_MPC = zmax*0;

% Implicitly assume that we apply U = 0 at the previous timestep since we
% start at the origin
v_i = controlArgs.v0;


% Set max iters
options.MaxIter = 10000;

for i=1:N-1
    % Get current state vector
    Xi    = X(:,i);

    % Execute control law at current time-step - If unspecified flag is
    % passed in, the dynamics will integrate with no control law
    ref_i = r(:,i);

    % Primal shifted warm-start - for compressed
    if i == 1
        xN = x0;
        kappa_i = 0;
        v_i = v0;
    else
        xPrev = X(:,i-1);
        xVec = AHat*xPrev + BHat*U_MPC;
        xN = xVec(end-n+1:end);
        uBar_i = Gu*v_i;
        xBar_i = Gx*v_i;
        uN = uBar_i - K*(xN - xBar_i);
        U_Shift = [U_MPC(m+1:end); uN];

        % Currrent psi
        xPert_i = Xi - Gx*v_i;
        c_QP = G_MPC*xPert_i;
        VPrime = U_Shift'*H_QP*U_Shift + 2*c_QP'*U_Shift + xPert_i'*W_MPC*xPert_i;
        psiPrime = sqrt(VPrime);

        % Select kappa
        if psiPrime > const.r_psi
            kappa_i = 0;
        else
            denomTerm = Gx*(ref_i - v_i);
            kappaPrime = (const.r_psi - psiPrime)/sqrt(denomTerm'*W_MPC*denomTerm);
            kappa_i = min([1, kappaPrime]);
        end
        v_i = v_i + kappa_i*(ref_i - v_i);
    end

    % Compute MPC solution
    xTilde_i = Xi - Gx*v_i;
    c_QP = G_MPC*xTilde_i;
    [U_MPC,iter_i] = projGradSolver(H_QP,c_QP,U_MPC*0,zmin,zmax,options);
    

    % Integrate dynamics
    Ui =  U_MPC(1:m);
    X(:,i+1) = A*Xi + B*Ui;
    U(:,i)   = Ui;

    % Store variables
    refHist(:,i) = v_i;
    cost_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
    V_i = cost_i;
    costHist(i) = cost_i;
    VHist(i) = V_i;
    errorHist(i) = 0;
    ellHist(i) = iter_i;
    ellBarHist(i) = -1;
    kappaHist(i) = kappa_i;
end

% Define outputs
output.X = X;
output.U = U;
output.refHist = refHist;
output.kappaHist = kappaHist;
output.VHist = VHist;
output.costHist = costHist;
output.errorHist = errorHist;
output.ellHist = ellHist;
output.ellBarHist = ellBarHist;
end % end function





