function [Optis, fun] = search(L,d,a,p)
%  
% [Optis,fun] = search(L,d,a,p) produces the p optimal solutions to the ILS
%             problem  min_{z}f(z)=(z-a)'*(L'DL)^{-1}*(z-a) by search.
%
% Inputs:
%    L - n by n unit lower triangular matrix
%    d - n-dimensional real positive vector, d = diag(D)
%    a - n-dimensional real vector
%    p - number of required optimal solutions with default value 1
%
% Outputs:
%    Optis - n by p matrix whose j-th column is the j-th optimal solution
%    fun -  p-dimensional real vector for the p-values of the objective 
%           function at the p optimal solutions, referred to as residuals 
%           or distances, i.e., fun(j)=f(Optis(:,j))

% ------------------------------------------------------------------------
% Main reference:
%   X.-W. Chang, X. Yang, and T. Zhou, MLAMBDA: A Modified LAMBDA Method 
%   for Integer Least-squares Estimation, Journal of Geodesy, 79 (2005), 
%   pp. 552-565. 
%
% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Taha Ghassemi, Xiaohu Xie
% Copyright (c) 2011-2016. Scientific Computing Lab, McGill University.
% September 2011. Last revision: April 2016
% ------------------------------------------------------------------------

n = length(d);

% Initialization 
z = zeros(n,1);        % current candidate solution
zb = zeros(n,1);       % used to compute z
step = zeros(n,1);     % step(k) is used to enumerate z(k) at level k
dist = zeros(n+1,1);   % partial distance of a candidate solution
S = zeros(n,n);        % used for computing zb
path = n*ones(n,1);    % path(k) used for updating S(k,k:path(k)-1) 
Optis = zeros(n,p);    % to store the p candidate solutions  
fun = zeros(1,p);      % residuals at the p candidate solutions
imax = p;
maxDist = inf;         % initial ellipsoid bound 
count = 0;             % initial number of candidate solutions
ulevel = 1;

% Determine which level to move to after z(1) is chosen at level 1.
if p == 1            
    k0 = 2;
else
    k0 = 1;
end

% Initialization at level n
zb(n) = a(n);
z(n) = round(zb(n));
y(n) = zb(n) - z(n);
step(n) = sign(y(n));
k = n;

endSearch = false;


while ~endSearch
    path(ulevel:k-1) = k;
    for j = ulevel-1 : -1 : 1
        if path(j) < k 
            path(j) = k;
        else
            break;
        end
    end
    
    newDist = dist(k) + y(k) * y(k) / d(k);
    
    while newDist < maxDist
        if k ~= 1      % move to level k-1
            k = k - 1;
            dist(k) = newDist;
            for j = path(k) : -1 : k + 1
                S(j-1,k) = S(j,k) - y(j) * L(j,k);
            end
            zb(k) = a(k) + S(k,k);
            z(k) = round(zb(k));
            y(k) = zb(k) - z(k);
            step(k) = sign(y(k));
        else      % store the found point and try next valid integer
            if count < p - 1     % store the first p-1 initial points
                count  =  count + 1;
                Optis(:, count) = z;
                fun(count) = newDist;  % store f(z)
            else
                Optis(:, imax) = z;
                fun(imax) = newDist;
                [maxDist, imax] = max(fun);
            end
            k = k0;
            z(k) = z(k) + step(k);   % next valid integer      
            y(k) = zb(k) - z(k);
            step(k) = -step(k) - sign(step(k));
        end
        newDist = dist(k) + y(k) * y(k) / d(k);
    end
    ulevel = k;
    
    while newDist >= maxDist   % exit or move to level k+1
        if k == n             
            endSearch = true;
            break;
        end
        k = k + 1;      % move to level k+1
        z(k) = z(k) + step(k);    % next integer
        y(k) = zb(k) - z(k);
        step(k) = -step(k) - sign(step(k));
        newDist = dist(k) + y(k) * y(k) / d(k);
    end
end

% Sort the solutions by their corresponding residuals
[fun,p] = sort(fun);
Optis = Optis(:,p);

