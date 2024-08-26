function [eigVal] = weightEig(M,P,flag)
% Computes the eigvalues of lambda^P(M)

sqrtP = sqrtm(P);
switch flag
    case '+'
        eigVal = max(eig(sqrtP\(M/sqrtP)));
    case '-'
        eigVal = min(eig(sqrtP\(M/sqrtP)));
    otherwise 
        errror('Input + or - for flag')
end

