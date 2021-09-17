function [data]=getSatellitesInfo(epochRover,epochBase,nav,p0,pbase,Na,sats,elevMask)
% Get all info about available sats in the k-th epoch
%% constants
wie=7292115.1467e-11; %Earth rotation rate (rad/s)
c=2.99792458e8; %Speed of light in vacuum (m/s)
f1=1.57542*1e9; %L1 GPS frequency
lambda1=c/f1; %L1 GPS wavelength
%%
time=epochRover(1,2); %toe
Ni=size(epochRover,1); %number of satellites (=measurements)
%% allocate memmory
flagMissSat=zeros(Ni,1); %sats that are not found in the base obs
psat=zeros(3,Ni);
elev=zeros(Ni,1);
azim=zeros(Ni,1);
dPhi=zeros(Ni,1);
dRho=zeros(Ni,1);
hat_dRho=zeros(Ni,1);
hat_dPhi=zeros(Ni,1);
hatN=zeros(Ni,1);
L=zeros(3,Ni); %LOS (ecef)
%%
lla=SingleLlaFromEcef(p0);
for m=1:Ni
    satId=epochRover(m,4); %satId
    baseObs=epochBase(epochBase(:,4)==satId,:);
    if size(baseObs,1)==1
        % sat found in the base obs
        id=find(nav(1,:)==satId);  %select all ephemeris for the satellite
        [~,j]=min(abs(time-nav(18,id)));  %seek for the most recent orbit parameters
        j=id(j);
        tems=time-epochRover(m,5)/c; %approx. emission time (GPST)
        psat(:,m)=satPosition(nav(:,j),tems);
        theta=wie*norm(psat(:,m)-p0)/c;
        psat(:,m)=rotZ(theta)*psat(:,m); %compensate Earth-rotation
        Cen=DCM_en(lla(1),lla(2));
        L(:,m)=(psat(:,m)-p0)/norm(psat(:,m)-p0);
        los=Cen'*L(:,m); %Line of Sight (LOS) (NED)
        azim(m)=atan2(los(1),los(2)); %Azimute (NED frame)
        elev(m)=asin(-los(3)); %Elevation (rad) (NED frame)
        %% Single-Difference Measurement
        dRho(m)=epochRover(m,5)-baseObs(5);
        dPhi(m)=lambda1*(epochRover(m,8)-baseObs(8));
        %% Ambiguity
        if epochRover(m,9)==0
            hatN(m)=Na(sats==satId);
        else
            hatN(m)=round((dPhi(m)-dRho(m))/lambda1);
        end
        %% Prediction
        hat_dRho(m)=norm(psat(:,m)-p0)-norm(psat(:,m)-pbase);
        hat_dPhi(m)=hat_dRho(m)+lambda1*hatN(m);
    else
        flagMissSat(m)=1;
    end
end
%% Sort
data=[epochRover(:,4) elev dRho dPhi flagMissSat psat' L' epochRover(:,9) hat_dRho hat_dPhi hatN];
data=sortrows(data,-2);
data(data(:,2)<deg2rad(elevMask),:)=[];
data(data(:,5)==1,:)=[];
%% result
% satInfo.satId=epochRover(:,4);
% satInfo.dRho=dRho;
% satInfo.dPhi=dPhi;
% satInfo.flagMissSat=flagMissSat;
% satInfo.elev=elev;
% satInfo.azim=azim;
% satInfo.psat=psat';
% satInfo.LLI=epochRover(:,9);
% satInfo.L=L';
% satInfo.hat_dRho=bar_dRho;
% satInfo.hat_dPhi=bar_dRho+lambda1*hatN;
% satInfo.hatN=hatN;
end