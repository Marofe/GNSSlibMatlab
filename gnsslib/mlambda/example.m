%
% A small example to run function mlambda.m
%

% Construct data
m = 5; n = 3; 
A = randn(m,n);
x_true = [1; -2; 3];
y = A*x_true + 1.e-3*randn(m,1);
a = A\y;
W = inv(A'*A)
p = 2;

% Find p optimal solutions to the ILS problem min_{x}(x-a)'*W^{-1}*(x-a)
[X, r] = mlambda(W,a,p);

display('The true integer parameter vector')
x_true

display('The two integer least squares estimates')
X

