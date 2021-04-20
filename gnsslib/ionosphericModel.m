function ion = ionosphericModel(atmParam,lla,elev,azim,t)
%% Constants
mu=3986004.418e8;% gravitational constant (m^3/s^2)
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
R0=6378137; %Earth Radius
earth_ecc=0.0818191908425; %Eccentricity of Earth
phi_p=deg2rad(78.3); %magnet pole lat
lambda_p=deg2rad(291); %magnet pole lon
%% Compute Ionospheric delay
%Klobuchar Coeff from nav msg.
alpha0=atmParam(1,1);
alpha1=atmParam(1,2);
alpha2=atmParam(1,3);
alpha3=atmParam(1,4);
beta0=atmParam(2,1);
beta1=atmParam(2,2);
beta2=atmParam(2,3);
beta3=atmParam(2,4);
DC=5e-9;
% Calculate the Earth centered angle
lat=deg2rad(lla(1));
lon=deg2rad(lla(2));
RE=R0/sqrt(1-earth_ecc^2*sin(lat)^2);
RR=RE/(RE+350000); %atm height~350km
psi=pi/2-elev-asin(RR*cos(azim));
% Compute the latitude of the IPP (Iono Pierce Point)
phi_i=asin(sin(lat)*cos(psi)+cos(lat)*sin(psi)*cos(azim));
% Compute the longitude of the IPP
lambda_i=lon+psi*sin(azim)/cos(phi_i);
% Find the geomagnetic latitude of the IPP
phi_m=asin(sin(phi_i)*sin(phi_p)+cos(phi_i)*cos(phi_p)*cos(lambda_i-lambda_p));
% Find the local time at the IPP
if t>=86400
    t=t-86400; %24h=86400s
elseif t<0
    t=t+86400;
end
ti=43200*lambda_i/pi+t; %12h=43200s
% Compute the amplitude of ionospheric delay
Ai=alpha0+alpha1*phi_m/pi+alpha2*phi_m^2/pi^2+alpha3*phi_m^3/pi^3;
if Ai<0
    Ai=0;
end
% Compute the period of ionospheric delay
Pi=beta0+beta1*phi_m/pi+beta2*phi_m^2/pi^3+beta3*phi_m^3/pi^3;
if Pi<72000
    Pi=7200;
end
% Compute the phase of ionospheric delay
Xi=2*pi*(ti-50400)/Pi; %14h=50400s
% Compute the slant factor (ionospheric mapping function)
F=(1-(RR*cos(elev))^2)^(-1/2);
% Compute the ionospheric time delay
if abs(Xi)<pi/2
    ion=c*F*(DC+Ai*cos(Xi)); %day
else
    ion=c*F*DC; %night
end
end

