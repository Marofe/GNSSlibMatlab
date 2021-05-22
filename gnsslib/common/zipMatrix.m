function m = zipMatrix(M,n)
%Compact any symmetric matrix to a vector of n*(n+1)/2 elements
U=triu(M);
m=zeros(1,n*(n+1)/2);
j=1;
for i=1:n
    m(j:j+i-1)=U(1:i,i);
    j=j+i;
end
end

