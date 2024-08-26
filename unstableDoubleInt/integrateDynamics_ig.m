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
const = controlArgs.const;
eqFlag = controlArgs.eqFlag;

% Generate the constant QP matrices (fixed horizon lengths)
[H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,S,M] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);
H_QP = H_MPC;
HInv = inv(H_MPC);
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
phiHat_i = 0;
kappa_i = 0;
v_i = controlArgs.v0;
xTilde_i = x0 - Gx*v_i;
Ui = zeros(m,1);
c_QP = G_MPC*xTilde_i;
lambdaRMax = max(eig(R));

% Initialize psiCheck
psiCheck_i = const.r_psi/const.gamma;
LCheck_i = (const.r_psi + const.zeta*const.r_phi)/const.epsilon;
phase1Flag = 1;

for i=1:N-1
    % Get current state vector
    Xi    = X(:,i);

    % Execute control law at current time-step - If unspecified flag is
    % passed in, the dynamics will integrate with no control law
    ref_i = r(:,i);

    % Run MPC
    if ws_flag && eqFlag
        % *************** WARM-STARTING CASE ***************
        if kappa_i < 1
            % Calc psiHat_{k-1}
            prevCost = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            psiHat_i = min([const.r_psi, sqrt(prevCost)]);

            % Calc psiHatPrime_{k}
            ellHat_i = xTilde_i'*Q*xTilde_i + max([0, Ui'*R*Ui - lambdaRMax*phiHat_i]);
            psiPrime_i = sqrt(psiHat_i^2 - ellHat_i) + const.mu*phiHat_i;

            % Calc psiBar_k
            xPert_i = Xi - Gx*v_i;
            c_QP = G_MPC*xPert_i;
            costSample_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xPert_i'*W_MPC*xPert_i;
            psiBar_i = min([psiPrime_i, sqrt(costSample_i), const.r_psi]);

            % Set kappa and v
            if psiBar_i <= const.r_psi*const.gamma
                denomTerm = Gx*(ref_i - v_i);
                kappaPrime = (const.r_psi - psiBar_i)/sqrt(denomTerm'*W_MPC*denomTerm);
                kappa_i = min([1, kappaPrime]);
                psiTilde_i = const.r_psi;
            else
                kappa_i = 0;
                psiTilde_i = psiBar_i;
            end

            % Set QP matrices
            v_i = v_i + kappa_i*(ref_i-v_i);
            xTilde_prev = xTilde_i;
            xTilde_i = Xi - Gx*v_i;
            DeltaXTilde_i = (xTilde_i-xTilde_prev);
            DeltaXTildeNorm_i = sqrt(DeltaXTilde_i'*G_MPC'*HInv*G_MPC*DeltaXTilde_i);
            c_QP = G_MPC*xTilde_i;

            % Calculate ellBar^1 and ellBar^2
            ellBar_i_1 =  log(const.gamma*const.r_psi - sqrt(psiTilde_i^2 - xTilde_i'*Q*xTilde_i))/log(const.eta) ...
                - log(const.mu*phiHat_i+ const.mu*const.b*DeltaXTildeNorm_i)/log(const.eta);
            ellBar_i_2 = (log(const.r_phi) - log(phiHat_i + const.b*DeltaXTildeNorm_i))/log(const.eta);
            ellBar_i = max([ellBar_i_1, ellBar_i_2]);
            if ellBar_i < const.ellStar_ws
                ell_i = ceil(ellBar_i);
            else
                ell_i = ceil(const.ellStar_ws);
            end

            % Set phiHat_i
            phiHat_i = min([const.r_phi, const.eta^(ell_i)*(phiHat_i + const.b*DeltaXTildeNorm_i)]);


            % Set max iters
            options.MaxIter = ell_i;
        else % kappa_i = 1 case, WARM-STARTING
            % Calc psiHat_{k-1}
            prevCost = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            if phase1Flag
                psiHat_i = min([const.r_psi, sqrt(prevCost)]);
                phase1Flag = 0;
            else
                psiHat_i = min([const.r_psi, sqrt(prevCost)]);
                % psiHat_i = min([psiBar_i, sqrt(prevCost)]);
            end

            % Calc psiHatPrime_{k}
            ellHat_i = xTilde_i'*Q*xTilde_i + max([0, Ui'*R*Ui - lambdaRMax*phiHat_i]);
            psiPrime_i = sqrt(psiHat_i^2 - ellHat_i) + const.mu*phiHat_i;

            % Calc psiBar_k and set QP matrices
            xTilde_prev = xTilde_i;
            xTilde_i = Xi - Gx*ref_i;
            DeltaXTilde_i = xTilde_i - xTilde_prev;
            DeltaXTildeNorm_i = sqrt(DeltaXTilde_i'*G_MPC'*HInv*G_MPC*DeltaXTilde_i);
            c_QP = G_MPC*xTilde_i;
            costSample_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            psiBar_i = min([psiPrime_i, sqrt(costSample_i), const.r_psi]);

            % Set LCheck and psiCheck
            LCheck_i = const.epsilon*LCheck_i;

            % Set iterations
            if psiBar_i < LCheck_i
                phiCheck_i = min([(LCheck_i - psiBar_i)/const.zeta, const.r_phi]);
                ellBar_i_1 =  log(const.r_psi - sqrt(psiBar_i^2 - xTilde_i'*Q*xTilde_i))/log(const.eta) ...
                    - log(const.mu*phiHat_i+ const.mu*const.b*DeltaXTildeNorm_i)/log(const.eta);
                ellBar_i_2 = (log(phiCheck_i) - log(phiHat_i + const.b*DeltaXTildeNorm_i))/log(const.eta);
                ellBar_i = max([ellBar_i_1, ellBar_i_2]);
                if ellBar_i < const.ellStar_ws
                    ell_i = ceil(ellBar_i);
                else
                    ell_i = ceil(const.ellStar_ws);
                end
            else
                ellBar_i = ceil(const.ellStar_ws);
                ell_i = ellBar_i;
            end

            % Set phiHat_i
            phiHat_i = min([const.r_phi, const.eta^(ell_i)*(phiHat_i + const.b*DeltaXTildeNorm_i)]);

            % Set max iters
            options.MaxIter = const.ellStar_ws;
        end

        % Compute MPC solution
        U_MPC = projGradSolver(H_QP,c_QP,U_MPC,zmin,zmax,options);

        % Compute the optimal solution
        options.MaxIter = 10000;
        U_OPT = projGradSolver(H_QP,c_QP,U_MPC,zmin,zmax,options);
        Error_i = norm(U_MPC - U_OPT);
    else % COLD-STARTING CASE
        % *************** COLD-STARTING CASE ***************
        if kappa_i < 1
            % Calc psiHat_{k-1}
            prevCost = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            psiHat_i = min([const.r_psi, sqrt(prevCost)]);

            % Calc psiHatPrime_{k}
            ellHat_i = xTilde_i'*Q*xTilde_i + max([0, Ui'*R*Ui - lambdaRMax*phiHat_i]);
            psiPrime_i = sqrt(psiHat_i^2 - ellHat_i) + const.mu*phiHat_i;

            % Calc psiBar_k
            xPert_i = Xi - Gx*v_i;
            c_QP = G_MPC*xPert_i;
            costSample_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xPert_i'*W_MPC*xPert_i;
            psiBar_i = min([psiPrime_i, sqrt(costSample_i), const.gamma*const.r_psi]);

            % Calc kappa
            denomTerm = Gx*(ref_i - v_i);
            kappaPrime = (const.r_psi - psiBar_i)/sqrt(denomTerm'*W_MPC*denomTerm);
            kappa_i = min([1, kappaPrime]);

            % Calc v_i
            v_i = v_i + kappa_i*(ref_i - v_i);

            % Set QP matrices
            xTilde_i = Xi - Gx*v_i;
            c_QP = G_MPC*xTilde_i;

            % Calc ellBar_k
            xTildeNorm_i = sqrt(xTilde_i'*G_MPC'*HInv*G_MPC*xTilde_i);
            ellBar_i = log(const.gamma*const.r_psi - sqrt(const.r_psi^2 - xTilde_i'*Q*xTilde_i))/log(const.eta) ...
                - log(const.mu*const.b*xTildeNorm_i)/log(const.eta);
            if ellBar_i < const.ellStar_cs
                ell_i = ceil(ellBar_i);
            else
                ell_i = ceil(const.ellStar_cs);
            end

            % Set phiHat_i
            phiHat_i = min([const.r_phi, const.eta^(ell_i)*const.b*xTildeNorm_i]);

            % Set max iters
            options.MaxIter = ell_i;

            % Switch eqFlag to 1 to get to the real WS case in the case
            % where wsFlag = 1
            eqFlag = 1;
        else
            % Calc psiHat_{k-1}, this is our estimate of the value
            % functuion at the PREVIOUS time step
            prevCost = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            if phase1Flag
                psiHat_i = min([const.r_psi, sqrt(prevCost)]);
                phase1Flag = 0;
                i_2 = 1;
            else
                psiHat_i = min([psiCheck_i, sqrt(prevCost)]);
                i_2 = i_2 + 1;
            end

            % Calc psiHatPrime_{k}
            ellHat_i = xTilde_i'*Q*xTilde_i + max([0, Ui'*R*Ui - lambdaRMax*phiHat_i]);
            psiPrime_i = sqrt(psiHat_i^2 - ellHat_i) + const.mu*phiHat_i;

            % Set QP matrices and set psiBar_i, psiBar_i is our estimate of
            % tbe value function at the CURRENT time step
            xTilde_i = Xi - Gx*ref_i;
            c_QP = G_MPC*xTilde_i;
            costSample_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
            psiBar_i = min([psiPrime_i,sqrt(costSample_i), const.gamma*psiHat_i]);

            % Set psiCheck_i, this is our TARGET value function for the
            % FOLLOWING time step
            psiCheck_i = const.gamma*psiCheck_i;

            % Calc ellBar_k
            xTildeNorm_i = sqrt(xTilde_i'*G_MPC'*HInv*G_MPC*xTilde_i);
            if psiCheck_i > sqrt(psiBar_i^2 - xTilde_i'*Q*xTilde_i)
                ellBar_i = log(psiCheck_i - sqrt(psiBar_i^2 - xTilde_i'*Q*xTilde_i))/log(const.eta) ...
                    - log(const.mu*const.b*xTildeNorm_i)/log(const.eta);
            else
                ellBar_i = 1;
            end
            if ~isreal(ellBar_i)
                ellBar_i
                error
            end
            if ellBar_i < const.ellStar_cs
                ell_i = max([ceil(ellBar_i),1]);
            else
                ell_i = ceil(const.ellStar_cs);
            end


            % Set phiHat_i
            phiHat_i = min([const.r_phi, const.eta^(ell_i)*const.b*xTildeNorm_i]);

            % Set max iters
            options.MaxIter = ell_i;
        end

        % Compute MPC solution
        U_MPC = projGradSolver(H_QP,c_QP,U_MPC*0,zmin,zmax,options);

        % Compute the optimal solution
        options.MaxIter = 100000;
        U_OPT = projGradSolver(H_QP,c_QP,U_MPC,zmin,zmax,options);
        Error_i = norm(U_MPC - U_OPT);
    end


    % Integrate dynamics
    Ui =  U_MPC(1:m);
    X(:,i+1) = A*Xi + B*Ui;
    U(:,i)   = Ui;

    % Store variables
    refHist(:,i) = v_i;
    cost_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
    V_i = U_OPT'*H_QP*U_OPT + 2*c_QP'*U_OPT + xTilde_i'*W_MPC*xTilde_i;
    costHist(i) = cost_i;
    VHist(i) = V_i;

    if psiBar_i^2 < V_i && kappa_i == 1 && phase1Flag == 0
        if (V_i - psiBar_i^2)/V_i < 0.01 
            psiBar_i = sqrt(V_i);
        else
             V_i
             psiBar_i
             i_2
        end
    end
    errorHist(i) = Error_i;
    ellHist(i) = ell_i;
    ellBarHist(i) = ellBar_i;
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


% output.terminalFlag = terminalFlag;
end % end function

% function [COST,terminalFlag] = evaluateCost(X0,U,N,m,A,B,P,Q,R,xBar,uBar,alpha,lambda)
%
% % Do the first point outside the loop to initialize
% X_i = X0;
% U_i = U(1:m);
% COST = (X_i - xBar)'*Q*(X_i-xBar) + (U_i - uBar)'*R*(U_i-uBar);
% X_i = A*X_i + B*U_i;
%
% % Do 2:N-1
% for i = 2:N
%     U_i = U(1 + (i-1)*m:m + (i-1)*m);
%     COST = COST + (X_i - xBar)'*Q*(X_i-xBar) + (U_i - uBar)'*R*(U_i-uBar);
%     X_i = A*X_i + B*U_i;
% end
%
% % Do N
% terminalCost = (X_i - xBar)'*P*(X_i-xBar);
% COST = COST + terminalCost;
%
% % Check terminal condition
% if nargin > 11
%     if terminalCost < lambda*alpha
%         terminalFlag = 1;
%     else
%         terminalFlag = 0;
%     end
% else
%     terminalFlag = [];
% end
%
% end




