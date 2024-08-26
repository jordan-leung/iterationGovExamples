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
% RUN_FLAG = 'ig';

% Set warm-start flag
% ws_flag = 1;
ws_flag = 0;

% Save figures and data?
saveFigFlag = 0;

% Run common file to set simulation parameters
setSimulationParameters;


%% Execute MPC solutions
NSample = 150;
x1Sample = linspace(-0.65,0.65,NSample); % 0 (to 1m) is threshold for stab at 0 angle
x2Sample = linspace(-0.2,0.2,NSample); % [-0.5,0.5 is threshold for stab at 0 angle
XMatrix_pos = zeros(n,NSample^2);
XMatrix_neg = zeros(n,NSample^2);
costPosVec = zeros(NSample^2,1);
costNegVec = zeros(NSample^2,1);
[XGrid,YGrid] = meshgrid(x1Sample,x2Sample);
costGrid = zeros(NSample,NSample);

% Set matrices
[H_MPC,G_MPC,W_MPC,ACon,FCon,LCon,AHat,BHat] = generateQPMatrices_compressed(N_MPC,A,B,lambda*P,Q,R,xmax,xmin,umax,umin);
H_QP = H_MPC;

% Generate minmax vector
zmax = zeros(size(H_MPC,1),1);
zmin = zmax;
for i = 1:N_MPC
    zmin(1 + m*(i-1) : m + m*(i-1)) = umin;
    zmax(1 + m*(i-1) : m + m*(i-1)) = umax;
end
options.MaxIter = 100000;
options.xTol = 1e-6;

% Loop
count_pos = 1;
count_neg = 1;
specialX = zeros(n,1);
specialCost = 1000;
for i1 = 1:NSample
    i1
    for i2 = 1:NSample
        % Get xTilde and c_QP
        xTilde_i = [x1Sample(i1); x2Sample(i2)];

        % Solve
        c_QP = G_MPC*xTilde_i;
        U_MPC = projGradSolver(H_QP,c_QP,zmin*0,zmin,zmax,options);

        % Calculate cost
        cost_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
        costGrid(i1,i2) = cost_i;

        % Get state vector
        xi_MPC = AHat*xTilde_i + BHat*U_MPC;

        % Get final state perturbation and final input
        xi_N = xi_MPC(end-(n-1):end);
        u_N = -K*xi_N;

        % Check the condition
        if u_N <= umax && u_N >= umin
            costPosVec(count_pos) = cost_i;
            XMatrix_pos(:,count_pos) = xTilde_i;
            count_pos = count_pos + 1;
        else
            costNegVec(count_neg) = cost_i;
            XMatrix_neg(:,count_neg) = xTilde_i;
            count_neg = count_neg+1;
            if cost_i < specialCost
                specialCost = cost_i;
                specialX = xTilde_i;
            end
        end
    end
end

% Trim
count_pos = count_pos - 1;
count_neg = count_neg - 1;
costPosVec = costPosVec(1:count_pos);
costNegVec = costNegVec(1:count_neg);
XMatrix_pos = XMatrix_pos(:,1:count_pos);
XMatrix_neg = XMatrix_neg(:,1:count_neg);

max(costPosVec)
min(costNegVec)
specialCost

%% Plot

% %{
% Plot positive points
figure
for i = 1:count_pos
    X_i = XMatrix_pos(:,i);
    plot(X_i(1),X_i(2),'g.','MarkerSize',1);
    if i == 1
        hold on; box on; grid on;
    end
end

% Plot negative points
for i = 1:count_neg
    X_i = XMatrix_neg(:,i);
    plot(X_i(1),X_i(2),'r.','MarkerSize',1);
end
% plot(specialX(1),specialX(2),'k.','markersize',15)

% Plot the sublevel set
VBar = 0.9*min(costNegVec);
figure(1)
hold on
contour(XGrid,YGrid,costGrid',[VBar VBar])
%}

% figure(2)
% surf(XGrid,YGrid,costGrid)
% view(0,90)

savestr = ['./Data/ROAData_N',num2str(N_MPC)];
save(savestr,'XMatrix_pos','XMatrix_neg','XGrid','YGrid','costGrid','costPosVec','costNegVec')




