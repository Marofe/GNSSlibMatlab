function [satInfo]=getSatellitesInfo(epochRover,epochBase,nav,p0)
% Get all info about available sats in the k-th epoch
%% constants
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
%%
time=epochRover(1,2); %toe
Ni=size(epochRover,1); %number of satellites (=measurements)
%% allocate memmory
flagMissSat=zeros(Ni,1); %sats that are not found in the base obs
psat=zeros(3,Ni);
elev=zeros(Ni,1);
azim=zeros(Ni,1);
L=zeros(3,Ni); %LOS (ecef)
%%
lla=SingleLlaFromEcef(p0);
for m=1:Ni
    satId=epochRover(m,4); %satId
    baseObs=epochBase(epochBase(:,4)==satId,:);
    if size(baseObs,1)==1
        % sat found in the base obs
        tems=time-epochRover(m,5)/c; %approx. emission time (GPST)
        id=find(nav(1,:)==satId);  %select all ephemeris for the satellite
        [~,j]=min(abs(time-nav(18,id)));  %seek for the most recent orbit parameters
        j=id(j);
        psat(:,m)=satPosition(nav(:,j),tems);
        theta=wie*norm(psat(:,m)-p0)/c;
        psat(:,m)=rotZ(theta)*psat(:,m); %compensate Earth-rotation
        Cen=DCM_en(lla(1),lla(2));
        L(:,m)=(psat(:,m)-p0)/norm(psat(:,m)-p0);
        los=Cen'*L(:,m); %Line of Sight (LOS) (NED)
        azim(m)=atan2(los(1),los(2)); %Azimute (NED frame)
        elev(m)=asin(-los(3)); %Elevation (rad) (NED frame)
    else
        flagMissSat(m)=1;
    end
end
%% result
satInfo=[epochRover(:,4) flagMissSat elev azim psat' L'];
end