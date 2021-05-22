function M = unzipMatrix(m,n)
% Convert vector of n*(n+1)/2 elements in a symmetric nxn matrix
U=zeros(n);
j=1;
for i=1:n
    U(1:i,i)=m(j:j+i-1)';
    j=j+i;
end
D=diag(diag(U));
M=U-D+U';
end

