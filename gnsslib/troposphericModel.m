function [trop]=troposphericModel(lla,D,elev)
param0=[1013.25 299.65 26.31 6.3e-3 2.77;...
    1017.25 294.15 21.79 6.05e-3 3.15;...
    1015.75 283.15 11.66 5.58e-3 2.57;...
    1011.75 272.15 6.78 5.39e-3 1.81;...
    1013.00 263.65 4.11 4.53e-3 1.55];
param=[0 0 0 0 0;...
    -3.75 7 8.85 0.25e-3 0.33;...
    -2.25 11 7.24 0.32e-3 0.46;...
    -1.75 15 5.36 0.81e-3 0.74;...
    -0.50 14.50 3.39 0.62e-3 0.30];
X=[15 30 45 60 75]';
lat=lla(1);
H=lla(3);
%% parameters
k1=77.604; %(K/mbar)
k2=382000; %(K^2/mbar)
Rd=287.054; %(J/kgK) specific constant for dry air
gm=9.784; %(m/s^2)
g=9.80665; %(m/s^2)

%% Interpolation
lat=abs(lat);
if lat>15 && lat<75
    Xi0=interp1(X,param0,lat)';
    dXi=interp1(X,param,lat)';
elseif lat<=15
    Xi0=param0(1,:)';
    dXi=param(1,:)';
else
    Xi0=param0(end,:)';
    dXi=param(end,:)';
end
%% Atm parameters
if lla(1)>0
    %Northern latitude
    Dmin=28;
else
    %Southern latitude
    Dmin=211;
end
Xi=Xi0-dXi*cos(2*pi*(D-Dmin)/365.25);
P=Xi(1);
T=Xi(2);
e=Xi(3);
beta=Xi(4);
lambda=Xi(5);
%% Zenith delay
Trzd0=1e-6*k1*Rd*P/gm;
Trzw0=(1e-6*k2*Rd/((lambda+1)*gm-beta*Rd))*e/T;
K=1-beta*H/T;
Trzd=(K^(g/(Rd*beta)))*Trzd0;
Trzw=(K^((lambda+1)*g/(Rd*beta)-1))*Trzw0;
%%
M=1.001/sqrt(0.002001+sin(elev)^2);
trop=(Trzd+Trzw)*M;
%% Simplified model
%         dwet=0.1;
%         ddry=2.3*exp(-0.116*1e-3*lla(3));
%         tilm=1.001/sqrt(0.002001+sin(elev)^2);
%         trop2=(ddry+dwet)*tilm; %~90% of the delay
end