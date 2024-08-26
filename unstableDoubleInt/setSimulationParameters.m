% * Time vector
dt = 0.1;
t = 0:dt:30;
N = length(t);

% * Initial Conditions
X0 = [0; 0];
n = length(X0);

% * Define commanded reference
nr = 1;

% Case 1 reference
rConst =  0;
r =  rConst*ones(nr,N);


% State Vector (X): x       = cart position (m)
%                   xDot    = cart speed (m/s)

Ac = [0, 1;...
    0, 0.5];

Bc = [0; 1];

% Get discrete state-space matrices
m = size(Bc,2);
[A,B] = c2d(Ac,Bc,dt);


%% Simulation and MPC parameters

% Min-max constraints
xmax = [Inf; 0.5]; % this dont rly do anything, just placeholder since we use PGM
xmin = -xmax;
umax = 0.1;
umin = -0.1;

% Get tracking matrices
Es = [1 0];
F = 0;
[nx,nu] = size(B);
Z = [eye(nx)-A, B, zeros(nx,nu); Es, F, -eye(nr)];
null_Z = null(Z);
Gx = null_Z(1:nx,:);
Gu = null_Z(nx+1:nx+nu,:);
Gr = null_Z(nx+nu+1:end,:);
if size(Gr,1) == size(Gr,2)
    if det(Gr) ~= 0
        Gx = Gx*inv(Gr);
        Gu = Gu*inv(Gr);
        Gr = eye(nr);
    end
end


% Set horizon
N_MPC = 80;
T_MPC = dt;
lambda = 1;

% * Set P, Q, R
Q = 1*eye(n);
Q(1,1) = 1;
R = 1; % Control effort cost


% % Set Q and R for the terminal controller
Q_LQR = Q;
R_LQR = 1;
K = dlqr(A,B,Q_LQR,R_LQR);
lambda_1 = 1;
P = dlyap((A-B*K)',lambda_1*(Q + K'*R*K));
P = (P+P')/2;

%% Calculate ROA size

% -------------- Calculate constant d and alpha --------------
% Compute alpha according to Marco's formula - taking the minimum
% OmegaBar = { x | -Kx in U }
ACon = [-K;...
    K];
bCon = [umax
    umax];
alpha = 1e12; % initialize
invP = inv(P);
for i = 1:length(bCon)
    b_i = bCon(i);
    a_i = ACon(i,:)';
    VNew = (-b_i)^2/(a_i'*invP*a_i); % assumes x_eq = 0
    if VNew < alpha
        alpha = VNew;
    end
end
sqrtQ  = sqrtm(Q);
invsqrtQ = inv(sqrtQ);
eigVec = eig(invsqrtQ*P*invsqrtQ);
dConst = alpha/max(eigVec);
VBar_theory = N_MPC*dConst + lambda*alpha;