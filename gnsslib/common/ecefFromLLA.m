function [x,y,z] = ecefFromLLA(lat,lon,alt)
%ECEFFROMLLA Summary of this function goes here
% input->lat,lon,alt in rad and meters
% output-> x,y,z in meters
lat=deg2rad(lat);
lon=deg2rad(lon);
R0=6.378137e6;
e=0.0818191908425;
RN=R0*(1-e^2)./((1-e^2*sin(lat).^2).^(1/2));
RE=R0./sqrt(1-e^2*sin(lat).^2);
x=(RE+alt).*cos(lat).*cos(lon);
y=(RE+alt).*cos(lat).*sin(lon);
z=(RN+alt).*sin(lat);
end

