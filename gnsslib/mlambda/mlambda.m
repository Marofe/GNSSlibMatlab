function [X,r] = mlambda(W,a,p)
%
% [X,r] = mlambda(W,a,p) produces p optimal solutions to the integer 
%         least squares problem min_{x}(x-a)'*W^{-1}*(x-a) and the 
%         corresponding values of the objective function 
%
% Inputs:
%    W - n by n symmetric positive definite matrix. In GNSS, it is  
%        the covariane matrix of the real least squares estimator  
%        of the integer ambiguity vector.
%    a - n-dimensional real vector. In GNSS, it is the real least 
%        squares estimator of the ambiguity vector.
%    p - number of required optimal solutions with default value 1.
%
% Outputs:
%    X - n by p matrix. Its j-th column is the j-th optimal solution,
%        i.e., its objective function is the j-th smallest.
%    r - p-dimensional real vector for the p values of the objective 
%        function at the p optimal solutions, i.e.,
%        r(j)=(X(:,j)-a)'*W^{-1}*(X(:,j)-a)
%

% ------------------------------------------------------------------------
% Main references:
% [1] M. Al Borno, X.-W. Chang, and X. Xie. On "Decorrelation" in Solving 
%     Integer Least-Squares Problems for Ambiguity Determination, 
%     Survey Review, 46:37-49, 2014. 
% [2] X.-W. Chang, X. Yang, and T. Zhou, MLAMBDA: A Modified LAMBDA Method 
%     for Integer Least-squares Estimation, Journal of Geodesy, 79 (2005), 
%     pp. 552-565. 
%
% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiaohu Xie
% Copyright (c) 2011-2016. Scientific Computing Lab, McGill University.
% September 2011. Last revision: April 2016
% ------------------------------------------------------------------------

% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!')
end

[m,n] = size(W);
if m ~= n || m ~= size(a,1) || size(a,2) ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

% Transform the ILS problem by reduction
[L,d,Z,a] = reduction(W,a);

% Search to find the p optimal solutions of the reduced ILS problem
[Optis,r] = search(L,d,a,p);

% Find the p optimal solutions of the original ILS problem by transformation
X = Z'*Optis;
