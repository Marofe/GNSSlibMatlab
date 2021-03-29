function lla = SingleLlaFromEcef(pe)
%ECEFFROMLLA Summary of this function
% input-> [x,y,z] in meters
% output->[lat,lon,alt] in degree and meters
% Reference:
% Karl Osen. Accurate Conversion of Earth-Fixed Earth-Centered Coordinates 
% to Geodetic Coordinates. [Research Report]
% Norwegian University of Science and Technology. 2017. 
x=pe(1);
y=pe(2);
z=pe(3);
lla=zeros(3,1);
a=6.378137e6; %R0
e=0.0818191908425;
l=e^2/2;
Hmin=e^12/4;
%% Enhanced Algorithm
w2=x.^2+y.^2;
m=w2/(a^2);
n=z.^2*(1-e^2)/(a^2);
p=(m+n-4*l^2)/6;
G=m.*n.*l^2;
H=2*p.^3+G;
assert(sum(H<Hmin)==0,'H<Hmin.. not feasible');
C=(H+G+2*sqrt(H.*G)).^(1/3)/(2^(1/3));
i=-(2*l^2+m+n)/2;
P=p.^2;
beta=i/3-C-P./C;
k=l^2*(l^2-m-n);
t=sqrt(sqrt(beta.^2-k)-(beta+i)/2)-sign(m-n).*sqrt(abs((beta-i)/2));
F=t.^4+2*i.*t.^2+2*l*(m-n).*t+k;
dF=4*t.^3+4*i.*t+2*l.*(m-n);
dt=-F./dF;
u=t+dt+l;
v=t+dt-l;
w=sqrt(w2);
lat=atan2(z.*u,w.*v);
dw=w.*(1-1./u);
dz=z.*(1-(1-e^2)./v);
lon=atan2(y,x);
lla(1)=rad2deg(lat);%latitude
lla(2)=rad2deg(lon);%longitude
lla(3)=sign(u-1).*sqrt(dw.^2+dz.^2);%altitude
end

