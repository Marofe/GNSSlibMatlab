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
file='gnss_logger_dat-2021-01-27-13-14-20';
ephemeris = readRinexNav([file '.nav']);
[p0,allObs]=readRinexObs([file '.obs']);
% orbit=readSP3file('igs21423.sp3');
%% Load RTKLIB files
diffGnss=loadGnssData(['data/' file '-diff.pos']);
singleGnss=loadGnssData(['data/' file '-single-gps.pos']);
singleGlabGnss=loadGlabGnssData('glab-2021-01-27.pos');
%%
nav=ephemeris.gpsEphemeris;
atmParam=ephemeris.ionosphericParameters;
satID=unique(nav(1,:));
gpsObs=allObs(:,3)==1;
obs=allObs(gpsObs,:);
elevMask=15;
%% Least Square Solution
tic
%[time,p,ns]=leastSquareGpsSolutionSP3(orbit,nav,obs,p0,elevMask,atmParam);
[time,pls,nsls,dopls,atmDelayls]=leastSquareGpsSolution(nav,obs,p0,elevMask,atmParam);
fprintf('\nElapsed time=%.3f min',toc/60);
%% Iterated-EKF
tic
[time,pkf,nskf,dopkf,atmDelaykf]=ekfGpsSolution(nav,obs,p0,elevMask,atmParam);
fprintf('\nElapsed time=%.3f min',toc/60);
%% 
prefdiff=diffGnss(1,2:4)'
err0=norm(prefdiff-p0)
errdiff=diffGnss(1:end-1,2:4)-pls(:,1:3);
err=singleGnss(1:end-1,2:4)-pls(:,1:3);
referr=diffGnss(:,2:4)-singleGnss(:,2:4);
norm(errdiff(1,:))
mse1=sqrt(mean(mean(errdiff.^2)));
mse2=sqrt(mean(mean(referr.^2)));
fprintf('\nRMSE: %.3f (m)\n',mse1)
fprintf('RMSE RTKLIB: %.3f (m)\n',mse2)
fprintf('Gain: %.2f%%\n',(mse2/mse1-1)*100)
figure
plot(errdiff)
title('Error w.r.t Diff')
%%
[lat, lon, alt]=llaFromEcef(pls(:,1),pls(:,2),pls(:,3));
[latkf, lonkf, altkf]=llaFromEcef(pkf(:,1),pkf(:,2),pkf(:,3));
figure
plot(diffGnss(:,7))
hold on
% plot(singleGnss(:,7))
plot(alt)
grid on
legend('RTKLIB-Diff','SPS-Matlab')
figure
plot(diffGnss(1:end-1,7)-alt)
hold on
plot(nskf,'*-')
plot(diffGnss(:,7)-singleGnss(:,7))
plot(singleGnss(:,9),'o-')
title('height error')
legend('Matlab','ns','RTKLIB','ns')
grid on
figure
subplot(2,1,1)
plot(diffGnss(1:end-1,7)-alt)
hold on
plot(2*dopls(:,3),'k--')
plot(-2*dopls(:,3),'k--')
grid on
legend('Height Error','+vdop','-vdop')
title('ILS')
subplot(2,1,1)
plot(diffGnss(1:end-1,7)-altkf)
hold on
plot(2*dopkf(:,3),'k--')
plot(-2*dopkf(:,3),'k--')
grid on
legend('Height Error','+vdop','-vdop')
title('EKF')
subtitle('Height Residue')
figure
plot(diffGnss(:,5),diffGnss(:,6),'-')
hold on
% plot(singleGnss(:,5),singleGnss(:,6),'x-')
plot(lat,lon,'-')
% plot(latkf,lonkf,'*-')
legend('RTKLIB-Diff','SPP-Matlab')
grid on
figure
plot(atmDelayls)
legend('ZTD','Ion')
grid on