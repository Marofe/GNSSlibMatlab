close all
clear all
clc
format long
addLibrary('gnsslib')
addpath('data')
%% To-do-list
% -> Glonass
% -> Carrier-phase positioning
% -> IMM-Kalman-filter
% -> DGPS
% -> Kinematic
% RMSE: 0.854 (m)
% RMSE RTKLIB: 1.599 (m)
% Gain: 87.28%
%% GPS Standard Precision Processing (NLS)
% file='circular';
% file='rect';
% file='helic';
% rover='data/gnss_logger_dat-2021-01-27-13-14-20';
% base='data/base-gnss_logger_dat-2021-01-27-09-34-28';
% ephemeris = readRinexNav([base '.nav']);
% [p0base,allObsBase]=readRinexObs([base '.obs']);
% [p0rover,allObsRover]=readRinexObs([rover '.obs']);
% %% Load RTKLIB files (ground-truth)
% diffGnss=loadGnssData([rover '-diff.pos']);
% singleGnss=loadGnssData([rover '-single-gps.pos']);
load('gnss-data-20210127.mat')
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
%% Least Square Solution (DGPS Carrier Float)
tic
[time2,pdiff2,nsdiff2,dopdiff2]=leastSquareDgpsCarrierSolution(nav,obsBase,obsRover,p0base,p0rover,elevMask,atmParam);
fprintf('\nElapsed time=%.3f min',toc/60);
%% 
err2=diffGnss(1:end-1,2:4)-pls(:,1:3); %Single
err3=diffGnss(1:end-1,2:4)-pdiff2(:,1:3); %DGPS-carrier
mse2=sqrt(mean(mean(err2.^2))); %Single
mse3=sqrt(mean(mean(err3.^2))); %DGPS-carrier
fprintf('RMSE-Single: %.3f (m)\n',mse2)
fprintf('RMSE-DGPS-carrier: %.3f (m)\n',mse3)
fprintf('Gain Carrier/Single: %.2f%%\n',(mse2/mse3-1)*100)
figure
plot(err3)
title('Error DGPS w.r.t RTKLIB-Diff')
%%
[lat2, lon2, alt2]=llaFromEcef(pls(:,1),pls(:,2),pls(:,3));
[lat3, lon3, alt3]=llaFromEcef(pdiff2(:,1),pdiff2(:,2),pdiff2(:,3));
%%
figure
plot(diffGnss(:,7))
hold on
plot(alt2,'-')
plot(alt3,'-')
grid on
legend('RTKLIB-Diff','Single','Kinematic')
%%
figure
plot(dNa')
title('Ambiguity')