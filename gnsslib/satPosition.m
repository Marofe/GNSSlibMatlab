function [satPos,E,a]=satPosition(nav,time)
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
we=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
%% Nav Message Parameters
toe=nav(18);    %time of ephemeris
M0=nav(3);      %mean anomaly
sqrta=nav(4);   %semi major-axis square-root
a=sqrta^2;        %semi-major axis
deltan=nav(5);  %mean motion difference
ecc=nav(6);     %satellite orbit eccentricity
omega=nav(7);   %argument of perigee
cuc=nav(8);     %cos lat arg correction
cus=nav(9);     %sin lat arg correction
crc=nav(10);    %cos orbital radius correction
crs=nav(11);    %sin orbitatl radius correction
i0=nav(12);     %inclination at reference epoch (toe)
doti0=nav(13);  %rate of inclination angle
cic=nav(14);    %cos inclination correction
cis=nav(15);    %sin inclination correction
Omega0=nav(16); %Ascending node's right ascension
dotOmega=nav(17); %Rate of node's right ascension
%% Elapsed time (<4h)
tk=time-toe;
if tk>302400
    tk=tk-604800;
elseif tk<-302400
    tk=tk+604800;
end

%% Compute the mean anomaly
M=M0+(sqrt(mu/a^3)+deltan)*tk;

%% Solve (iteratively) the Kepler equation
E=M;
iter=1;
iterMax=100;
while abs(E-M-ecc*sin(E))>1e-15 && iter<=iterMax
    %Newton's method
    iter=iter+1;
    if iter>iterMax
        error('Max Iteration exceeded!')
    end
    E=E-(E-ecc*sin(E)-M)/(1-ecc*cos(E));
end

%% Compute the true anomaly
v=atan2(sqrt(1-ecc^2)*sin(E),cos(E)-ecc);
phi=omega+v;

%% Compute the argument of latitude
uk=omega+v+cuc*cos(2*phi)+cus*sin(2*phi);

%% Compute the radial distance
rk=a*(1-ecc*cos(E))+crc*cos(2*phi)+crs*sin(2*phi);

%% Compute the inclination
ik=i0+doti0*tk+cic*cos(2*phi)+cis*sin(2*phi);

%% Compute the longitude of the ascending node
lambdak=Omega0+(dotOmega-we)*tk-we*toe;

%% Compute the coordinates in the TRS (ECEF) frame
C=rotZ(-lambdak)*rotX(-ik)*rotZ(-uk);
satPos=C*[rk;0;0];
end