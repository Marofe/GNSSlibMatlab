function [X,Wm,Wc] = sigmaPoints(x0,P0,alpha,beta,kappa)
if size(x0,1)==1
    x0=x0';
end
n=size(x0,1);
P0=0.5*(P0+P0'); %re-inforce symmetry
S=chol(P0)';
lambda=alpha*(n+kappa)-n;
t=sqrt(lambda+size(P0,1));
X=[x0 x0+t*S x0-t*S];
Wm=[lambda/(n+lambda) 1/(2*(n+lambda))*ones(1,2*n)];
Wc=[lambda/(n+lambda)+(1-alpha^2+beta) 1/(2*(n+lambda))*ones(1,2*n)];
end

