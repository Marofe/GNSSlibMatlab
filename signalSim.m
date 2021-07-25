close all
clear all
clc
%% GPS Signal Simulation
f0=10.23*1e6; %sat clock frequency
f1=154*f0;%L1=1575.42*1e6;
f2=120*f0;%L2=1227.60*1e6
fdata=50; %data msg freq
fca=f0/10; %coarse aquisition freq
c=299792458;
lambda1=c/f1;
lambda2=c/f2;
%% Carrier
dt=1/(50*f0);
t=0:dt:500*dt;
carr=sin(2*pi*f0*t);
% figure
% plot(t,carr,'g','linewidth',2)
%% PRBS
dtca=1/f0;
tca=0:dtca:50*dtca;
ca=2*prbs(4,51)-1;
figure
stairs(tca,ca,'linewidth',2)
hold on
T=0:1/(100*f0):1/f0;
carr_ca=sin(2*pi*f0*T);
s=ca'*carr_ca;
s=s';
tt=0:1/(100*f0):(numel(ca)-1)/f0;
s=s(1:length(tt));
figure
% stairs(tca,ca,'linewidth',2)
hold on
plot(tt,s,'r','linewidth',2)
figure
carr=repmat(carr_ca,[1,50]);
plot(tt,carr(1:length(tt)),'g','linewidth',2)
%% data
dtca=1/f0;
tca=0:dtca:30*dtca;
ca=prbs(4,31);
figure
stairs(tca,ca,'linewidth',2)
%%
theta=-6*pi:0.01:6*pi;
plot(theta,atan2d(1,theta))