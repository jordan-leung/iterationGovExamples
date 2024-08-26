function  [x,iterCount,lastRes,execTime] = projGradSolver(H,f,x0,xl,xu,opts)
% Solves the quadratic programming problem:
% min 0.5*x'*H*x + f'*x   subject to:  xl <= x <= xu
% using the projected gradient method...

% Check to make sure provided initial condition matches the size of the QP
% matrices
if length(x0) ~= size(H,2)
    error('Mismatch in the size of x0. Check provided x0')
end

% Unpack options
if isfield(opts,'MaxIter')
    MaxIter = opts.MaxIter;
else
    MaxIter = 1e5;
end
if isfield(opts,'xTol')
    xTol = opts.xTol;
else
    xTol = 1e-6;
end
if isfield(opts,'printFlag')
    printFlag = opts.printFlag;
else
    printFlag = 0;
end
if isfield(opts,'eigH')
    eigH = opts.eigH;
else
    eigH = eig(H);
end


% Initialize
x = x0; % Initialize x
n = length(x);
iterCount = 0;



% Find the Lipshitz constant L and set step-size
L = max(eigH);
ell = min(eigH);
t = 2/(L + ell);


outerLoopRunCond = 1;
while outerLoopRunCond
    iterCount = iterCount + 1;

    % Store previous x value to check the convergence tolerance
    xPrev = x;
    
    % Current gradient direction
    g = H*x + f;
    
    % Compute temporary step
    xStar = x - t*g;
    
    % Project any components exceeding the constraints back onto the box
    % constraint
    for i = 1:n
        if xStar(i) > xu(i)
            x(i) = xu(i);
        elseif xStar(i) < xl(i)
            x(i) = xl(i);
        else
            x(i) = xStar(i);
        end
    end
    
    % Check convergence criteria
    lastRes = norm(x - xPrev,2);
    xConv = lastRes;
    
    if iterCount >= MaxIter
%         if printFlag
%             fprintf('-------------------------------------------------- \n')
%             fprintf(' *** Maximum iteration count of %0.0i exceeded. PG Algorithm Complete *** \n',MaxIter)
%             fprintf(' *** Max update on completition is dx = %0.2e *** \n',xConv)
%             fprintf('-------------------------------------------------- \n')
%         end
        warning('PGM exceeded maximum number of iterations')
        outerLoopRunCond = 0;
    elseif xConv < xTol
%         fprintf('-------------------------------------------------- \n')
%         fprintf(' *** max(dx) < tolx = %0.2e. PG Algorithm Complete *** \n',xTol)
%         fprintf(' *** Number of iteration for completition = %0.0i *** \n',iterCount)
%         fprintf('-------------------------------------------------- \n')
        outerLoopRunCond = 0;
    end
    
    if printFlag    
        fprintf('* Iteration: %0.0i, Residual: %0.02e \n',iterCount,xConv)
    end
end


end

