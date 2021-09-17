function [hx,P] = filterPred(hx0,P0,dt,Qve)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A=[eye(3) dt*eye(3);zeros(3) eye(3)];
hx=A*hx0;
Q=[Qve*dt^3/3 Qve*dt^2/2;Qve*dt^2/2 Qve*dt];
P=A*P0*A'+Q;
end

