dt = 0.025;
% t = 0:dt:15;
t = 0:dt:10;
N = length(t);

% * Define commanded reference
nr = 1;

% Case 1 reference
rConst =  0;
r =  rConst*ones(nr,N);


% State Vector (X): x       = cart position (m) 
%                   phi     = pendulum angle pertubation (theta - pi) (rad)
%                   xDot    = cart speed (m/s)
%                   phiDot  = pendulum angular velocity (rad/s)

% Continuous-time dynamics
% mCart  = 1;      % kg, mass of cart
mCart = 0.5;
% mPend  = 1;      % kg, mass of pendulum
mPend  = 0.1;      % kg, mass of pendulum
bf     = 1;      % N/m/s, damping/fraction coefficient --- nominally 1
ell    = 1;      % m, length to pendulum centre of mass from bolt point
Im     = 1/3*mPend*ell^2;    % kg*m^2, mass moment of inertia of pendulum
g      = 9.81;     % m/s^2, gravity
denom = Im*(mPend + mCart) + mPend*mCart*ell^2;
Ac = [0,                          1,                                 0,   0;...
     0, -(Im+mPend*ell^2)*bf/denom,             mPend^2*g*ell^2/denom,   0;...
     0,                          0,                                 0,   1;...
     0,        -mPend*ell*bf/denom,   mPend*g*ell*(mPend+mCart)/denom,   0];

Bc = [                     0;...
     (Im+mPend*ell^2)/denom;...
                          0;...
            mPend*ell/denom];
        

% * Initial Conditions
X0 = [2; 0; 0; 0]; % this is overwritten in the scripts where this is called
n = length(X0);

        
% Get discrete state-space matrices
m = size(Bc,2);
[A,B] = c2d(Ac,Bc,dt);


%% Simulation and MPC parameters

% Min-max constraints
xmax = 1000*ones(n,1); % this dont rly do anything, just placeholder since we use PGM
xmin = -xmax;
umax = 1;
umin = -umax;

% Get tracking matrices
Es = [1 0 0 0];
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
N_MPC = 35;
T_MPC = dt;
lambda = 1;

% * Set P, Q, R
Q = 1*eye(n);
Q(1,1) = 10;
Q(2,2) = 1;
Q(3,3) = 1;
Q(4,4) = 1;
R = 0.1; % Control effort cost

% % * Set P, Q, R
% Q = 1*eye(n);
% R = 0.1; % Control effort cost

% % Set Q and R for the terminal controller
% Q_LQR = eye(n);
% R_LQR = 1;
Q_LQR = Q;
R_LQR = R;
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
