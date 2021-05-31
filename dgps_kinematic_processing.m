close all
clear all
clc
format long
addLibrary('gnsslib')
addpath('data')
%% To-do-list
% -> Glonass
% -> IMM-Kalman-filter
% -> Kinematic: 
%           - add filter const. velocity
%           - pick as sat ref the highest elevation
%           - solve ambiguity
%           - smoothing
% RMSE: 0.854 (m)
% RMSE RTKLIB: 1.599 (m)
% Gain: 87.28%
%% GPS Standard Precision Processing (NLS)
% file='circular';
% file='rect';
% file='helic';
rover='data/gnss_logger_dat-2021-01-27-13-14-20';
base='data/base-gnss_logger_dat-2021-01-27-09-34-28';
ephemeris = readRinexNav([base '.nav']);
[p0base,allObsBase]=readRinexObs([base '.obs']);
[p0rover,allObsRover]=readRinexObs([rover '.obs']);
%% Load RTKLIB files (ground-truth)
diffGnss=loadGnssData([rover '-diff.pos']);
singleGnss=loadGnssData([rover '-single-gps.pos']);
%% Select only GPS data
nav=ephemeris.gpsEphemeris;
atmParam=ephemeris.ionosphericParameters;
satID=unique(nav(1,:));
gpsObsBase=allObsBase(:,3)==1;
obsBase=allObsBase(gpsObsBase,:);
gpsObsRover=allObsRover(:,3)==1;
obsRover=allObsRover(gpsObsRover,:);
%% Elevation Mask
elevMask=15;
%% Base coordinate
antHeight=2.105;
lat=[-23 17 37.9762];
lon=[-48 35 40.5226];
alt=623+antHeight;
[x,y,z]=ecefFromLLA(dms2degrees(lat),dms2degrees(lon),alt);
p0base=[x y z]';
%% Sky Plot
skyPlot(obsRover,nav,p0rover,elevMask)
%% Least Square Solution (Single)
% tic
% [time0,pls,nsls,dopls]=leastSquareGpsSolution(nav,obsRover,p0rover,elevMask,atmParam);
% fprintf('\nElapsed time=%.3f min',toc/60);
% save single time0 pls nsls dopls
load('single')
%% Least Square Solution (DGPS)
% tic
% [time1,pdiff,nsdiff,dopdiff]=leastSquareDgpsSolution(nav,obsBase,obsRover,p0base,p0rover,elevMask,atmParam);
% fprintf('\nElapsed time=%.3f min',toc/60);
% save dgps time1 pdiff nsdiff dopdiff
load('dgps')
%% Least Square Solution (Kinematic Float)
tic
[time2,pdiff2,nsdiff2,dopdiff2,dNa]=leastSquareKinematicFloatSolution(nav,obsBase,obsRover,p0base,p0rover,elevMask,atmParam);
fprintf('\nElapsed time=%.3f min',toc/60);
%% 
err1=diffGnss(1:end-1,2:4)-pdiff(:,1:3); %DGPS
err2=diffGnss(1:end-1,2:4)-pls(:,1:3); %Single
err3=diffGnss(1:end-1,2:4)-pdiff2(:,1:3); %DGPS-carrier
referr=diffGnss(:,2:4)-singleGnss(:,2:4); %Single-RTKLIB
mse1=sqrt(mean(mean(err1.^2))); %DGPS
mse2=sqrt(mean(mean(err2.^2))); %Single
mse3=sqrt(mean(mean(err3.^2))); %DGPS-carrier
mse0=sqrt(mean(mean(referr.^2))); %Single-RTKLIB
fprintf('\nRMSE-DGPS: %.3f (m)\n',mse1)
fprintf('RMSE-Single: %.3f (m)\n',mse2)
fprintf('RMSE-DGPS-carrier: %.3f (m)\n',mse3)
fprintf('RMSE Single-RTKLIB: %.3f (m)\n',mse0)
fprintf('Gain DGPS/Single: %.2f%%\n',(mse2/mse1-1)*100)
fprintf('Gain Carrier/DGPS: %.2f%%\n',(mse1/mse3-1)*100)
fprintf('Gain DGPS/RTKLIB: %.2f%%\n',(mse2/mse0-1)*100)
figure
plot(err1)
title('Error DGPS w.r.t RTKLIB-Diff')
%%
[lat1, lon1, alt1]=llaFromEcef(pdiff(:,1),pdiff(:,2),pdiff(:,3));
[lat2, lon2, alt2]=llaFromEcef(pls(:,1),pls(:,2),pls(:,3));
[lat3, lon3, alt3]=llaFromEcef(pdiff2(:,1),pdiff2(:,2),pdiff2(:,3));
figure
plot(diffGnss(:,7))
hold on
plot(alt1)
plot(alt2)
plot(alt3,'--')
grid on
legend('RTKLIB-Diff','DGPS','Single','DGPS-Carrier')
% figure
% plot(diffGnss(:,9),'>-')
% hold on
% plot(singleGnss(:,9),'<-')
% plot(nsdiff,'o-')
% plot(nsls,'s-')
% title('Number of Satellites')
% legend('RTKLIB-diff(GNSS)','RTKBLI-single(GPS)','DGPS','Single')
% grid on
figure
plot(diffGnss(1:end-1,7)-alt1)
hold on
plot(diffGnss(1:end-1,7)-alt2)
plot(diffGnss(1:end-1,7)-alt3)
plot(diffGnss(:,7)-singleGnss(:,7))
plot(2*dopls(:,3),'k--')
plot(-2*dopls(:,3),'k--')
grid on
legend('DGPS','Single','Single-RTKLIB','DGPS-carrier','+vdop','-vdop')
title('Height Error')
figure
plot(diffGnss(:,5),diffGnss(:,6),'o-')
hold on
plot(lat3,lon3,'*-')
legend('RTKLIB-Diff','DGPS-carrier')
grid on