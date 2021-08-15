close all
clear all
clc
format long
addLibrary('C:\Users\Marcos\Dropbox\Matlab\GNSS\gnsslib')
addpath('data')
%% To-do-list
% -> Glonass
% -> IMM-Kalman-filter
% -> UKF
% -> Kinematic: 
%           - improve the ambiguity correlation
%           - backward
%           - smoothing
%% GPS Standard Precision Processing (NLS)
% file='circular';
% file='rect';
% file='helic';
rover='data/gnss_logger_dat-2021-01-27-13-14-20';
% base='data/base-gnss_logger_dat-2021-01-27-09-34-28';
% ephemeris = readRinexNav([base '.nav']);
% [~,allObsBase]=readRinexObs([base '.obs']);
% [p0rover,allObsRover]=readRinexObs([rover '.obs']);
% %% Load RTKLIB files (ground-truth)
%diffGnss=loadGnssData([rover '-diff.pos']);
diffGnssFwd=loadGnssData('gnss_logger_dat-2021-01-27-13-14-20-sol1_fix_fwd.pos');
diffGnssBck=loadGnssData('gnss_logger_dat-2021-01-27-13-14-20-sol2_fix_bck.pos');
load('gnss-data-20210127.mat')
diffGnss=loadGnssData('gnss_logger_dat-2021-01-27-13-14-20-sol3_fix_fwd_bck_gps.pos');
%% Select GPS only data
gpsObsBase=allObsBase(:,3)==1;
obsBase=allObsBase(gpsObsBase,:);
gpsObsRover=allObsRover(:,3)==1;
obsRover=allObsRover(gpsObsRover,:);
%% Elevation Mask
elevMask=12;
%% Base coordinate
antHeight=2.105;
lat=[-23 17 37.9762];
lon=[-48 35 40.5226];
alt=623+antHeight;
[x,y,z]=ecefFromLLA(dms2degrees(lat),dms2degrees(lon),alt);
p0base=[x y z]';
%% Code vs Phase
c=2.99792458e8; %Speed of light in vacuum (m/s)
f1=1.57542*1e9; %L1 GPS frequency
lambda1=c/f1; %L1 GPS wavelength
%% clean obs data
obsRover(obsRover(:,8)==0,:)=[];
obsBase(obsBase(:,8)==0,:)=[];
obs=obsRover(obsRover(:,4)==obsRover(1,4),:);
obsRef=obsBase(obsBase(:,4)==obsRover(1,4),:);
%%
% figure
% plot(obs(:,2),obs(:,5))
% hold on
% plot(obsRef(:,2),obsRef(:,5))
% legend('Rover','Base')
% % plot(obs(:,2),obs(:,8)*lambda1,'--')
% grid on
%%
% figure
% plot(obs(:,2),obs(:,5)-obs(:,8)*lambda1)
% hold on
% obs(obs(:,9)==0,9)=NaN;
% obs(obs(:,6)==0,6)=NaN;
% stem(obs(:,2),obs(:,9),'r') %phase LLI
% stem(obs(:,2),obs(:,6),'b') %code LLI
%% Sky Plot
% skyPlot(obsRover,nav,p0rover,elevMask)
%% Least Square Solution (Single)
% tic
% [time0,pls,nsls,dopls]=leastSquareGpsSolution(nav,obsRover,p0rover,elevMask,atmParam);
% fprintf('\nElapsed time=%.3f min',toc/60);
% save single time0 pls nsls dopls
load('single')
%% Least Square Solution (Kinematic)
tic
[time2,pdiff2,dNa,satsOnView,sats,Ra,cycleSlip,res]=KinematicSolution(ephemeris,obsBase,obsRover,p0base,p0rover,elevMask);
dNa(dNa==0)=NaN;
res.rho(res.rho==0)=NaN;
res.phi(res.phi==0)=NaN;
fprintf('\nElapsed time=%.3f min',toc/60);
%% 
err2=diffGnss(1:end-1,2:4)-pls(:,1:3); %Single
err3=diffGnss(1:end-1,2:4)-pdiff2(:,1:3); %DGPS-carrier
mse2=sqrt(mean(mean(err2.^2))); %Single
mse3=sqrt(mean(mean(err3.^2))); %DGPS-carrier
fprintf('\nRMSE-Single: %.3f (m)\n',mse2)
fprintf('RMSE-DGPS-carrier: %.3f (m)\n',mse3)
fprintf('Gain Carrier/Single: %.2f%%\n',(mse2/mse3-1)*100)
figure
plot(err3)
title('Error DGPS w.r.t RTKLIB-Diff')
%%
[lat2, lon2, alt2]=llaFromEcef(pls(:,1),pls(:,2),pls(:,3));
[lat3, lon3, alt3]=llaFromEcef(pdiff2(:,1),pdiff2(:,2),pdiff2(:,3));
% [lat3s, lon3s, alt3s]=llaFromEcef(pdiffs(:,1),pdiffs(:,2),pdiffs(:,3));
%%
figure
plot(diffGnss(:,7),'linewidth',1.5)
hold on
plot(diffGnssFwd(:,7),'linewidth',1.5)
% plot(alt2,'-','linewidth',1.5)
plot(alt3,'-','linewidth',1.5)
plot(diffGnss(1,7)+(Ra>3),'s-g')
grid on
plot(diffGnss(1,7)+10*cycleSlip,'r-s')
% legend('RTKLIB-Diff','RTKLIB-fwd','Kinematic','FIX','Cycle-Slip')
%%
figure
plot(dNa')
hold on
stem(100*cycleSlip,'r')
title('Ambiguity')
%% Sats On View
% satsOnView(satsOnView==0)=NaN;
% figure
% hold on
% for j=1:numel(sats)
%     plot(satsOnView(j,:)*sats(j),'linewidth',2)
% end
figure
plot(Ra)
hold on
plot(sum(satsOnView))
%% Residues
figure,plot(res.rho'),title('CODE Residue')
figure,plot(res.phi'),title('PHASE Residue'),hold on,stem(cycleSlip,'r')