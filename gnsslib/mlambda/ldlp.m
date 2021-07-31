function [L,d,P,a] = ldlp(W,a)
%
% [L,d,P,a] = ldlp(W,a) computes the L'DL factorization of W with minimum
%             symmetric piovting: W = P'L'DLP and computes a:=Pa.
%
% Inputs:
%    W - n by n symmetric positive definite matrix. 
%    a - n-dimensional vector.
%
% Outputs:
%    L - n by n unit lower triangular matrix.
%    d - vector formed by the diagonal of the n by n diagonal matrix D.
%    P - n by n permutation matrix.
%    a - updated a.

% ------------------------------------------------------------------------
% Main references:
% [1] X.-W. Chang, X. Yang, and T. Zhou, MLAMBDA: A Modified LAMBDA Method 
%     for Integer Least-squares Estimation, Journal of Geodesy, 79 (2005), 
%     pp. 552-565. 
%
% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiaohu Xie 
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% September 2011. Last revision: April 2016
% ------------------------------------------------------------------------


n = size(W,1);
perm = 1:n;  % record the permutation information

for k = n:-1:1
    [mindiag, q] = min(diag(W(1:k,1:k)));
	if mindiag <= 0
		error ('Input matrix W is not positive definite!');
    end
    perm([k q]) = perm([q k]);
    W([k q],:) = W([q k],:);
    W(:,[k q]) = W(:, [q k]);
    W(k,1:k-1) = W(k,1:k-1)/W(k,k);
    W(1:k-1,1:k-1) = W(1:k-1,1:k-1) - W(k,1:k-1)'*W(k,k)*W(k,1:k-1);   
    a([k q]) = a([q k]);
end

d = diag(W)';
L = tril(W,-1) + eye(n);
P = zeros(n);
for j = 1:n
    P(j,perm(j)) = 1;
end

