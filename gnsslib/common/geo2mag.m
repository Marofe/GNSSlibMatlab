function mphi=geo2mag(lat,lon,elev,azim)
phi_pole=deg2rad(79); %latitude of geomagnetic pole
lam_pole=deg2rad(-71);% longitude of geomagnetic pole
lon=deg2rad(lat);
lat=deg2rad(lon);
Psi=0.0137/(elev+0.11)-0.022;
iphi=lat+Psi*cos(azim);
ilambda=lon+Psi*sin(azim)/cos(iphi);
mphi=iphi+0.064*cos(ilambda-1.617);
%mphi=rad2deg(acos(sin(lat)*sin(phi_pole)+cos(lat)*cos(phi_pole)*cos(lon-lam_pole)));
end