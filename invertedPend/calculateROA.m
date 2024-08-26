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
ws_flag = 1;
% ws_flag = 0;

% Save figures and data?
saveFigFlag = 0;

% Run common file to set simulation parameters
setSimulationParameters;


%% Execute MPC solutions
NSample = 15;
% NSample = 25;
x1Sample = linspace(0,1,NSample); % 0 (to 1m) is threshold for stab at 0 angle
x2Sample = linspace(-1.5,1.5,NSample); % [-0.5,0.5 is threshold for stab at 0 angle
x3Sample = linspace(-7.5*pi/180,7.5*pi/180,NSample); % 0 (to 1m) is threshold for stab at 0 angle
x4Sample = linspace(-0.2,0.2,NSample); % [-0.5,0.5 is threshold for stab at 0 angle
costPosVec = zeros(NSample^4,1);
costNegVec = zeros(NSample^4,1);
[XGrid,YGrid] = meshgrid(x1Sample,x3Sample);
costGrid = zeros(NSample,NSample);
feasGrid = zeros(NSample,NSample);

iMid = ceil(NSample/2);
if x4Sample(iMid)~=0
    error('not right mid')
end

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
options.MaxIter = 1000000;
options.xTol = 1e-6;

% Loop
count_pos = 1;
count_neg = 1;
for i1 = 1:NSample
    i1
    for i2 = 1:NSample
        for i3 = 1:NSample
            for i4 = 1:NSample
                % Get xTilde and c_QP
                xTilde_i = [x1Sample(i1); x2Sample(i2); x3Sample(i3); x4Sample(i4)];

                % Solve
                c_QP = G_MPC*xTilde_i;
                U_MPC = projGradSolver(H_QP,c_QP,zmin*0,zmin,zmax,options);

                % Calculate cost
                cost_i = U_MPC'*H_QP*U_MPC + 2*c_QP'*U_MPC + xTilde_i'*W_MPC*xTilde_i;
                if i2 == iMid && i4 == iMid
                    costGrid(i1,i3) = cost_i;
                end

                % Get state vector
                xi_MPC = AHat*xTilde_i + BHat*U_MPC;

                % Get final state perturbation and final input
                xi_N = xi_MPC(end-(n-1):end);
                u_N = -K*xi_N;

                % Check the condition
                if u_N <= umax && u_N >= umin
                    costPosVec(count_pos) = cost_i;
                    if i2 == iMid && i4 == iMid
                        feasGrid(i1,i3) = 1;
                    end
                    count_pos = count_pos + 1;
                else
                    costNegVec(count_neg) = cost_i;
                    count_neg = count_neg+1;
                end
            end
        end
    end
end

% Trim
count_pos = count_pos - 1;
count_neg = count_neg - 1;
costPosVec = costPosVec(1:count_pos);
costNegVec = costNegVec(1:count_neg);

max(costPosVec)
min(costNegVec)

%% Plot

% Plot positive points
figure
hold on; box on; grid on;
for i1 = 1:NSample
    for i3=1:NSample
        if feasGrid(i1,i3)
            plot(x1Sample(i1),x3Sample(i3),'g.','MarkerSize',1);
        else
            plot(x1Sample(i1),x3Sample(i3),'r.','MarkerSize',1);
        end
    end
end


% Plot the sublevel set
VBar = 0.8*min(costNegVec);
contour(XGrid,YGrid,costGrid',[VBar VBar])

% figure(2)
% surf(XGrid,YGrid,costGrid)
% view(0,90)

savestr = ['./Data/ROAData_N',num2str(N_MPC)];
% savestr = ['./Data/ROAData_N',num2str(N_MPC),'_final'];
save(savestr,'x1Sample','x2Sample','x3Sample','x4Sample','XGrid','YGrid','costGrid','feasGrid','costPosVec','costNegVec')

% Choose a maximum x2 appears to be 0.075



