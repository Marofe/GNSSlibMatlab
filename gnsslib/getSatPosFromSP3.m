function [pos,dt]=getSatPosFromSP3(orbit,time)
[~,j]=min(abs(orbit(:,1)-time));
x=interpLagran(orbit(j-5:j+5,1),orbit(j-5:j+5,2),time);
y=interpLagran(orbit(j-5:j+5,1),orbit(j-5:j+5,3),time);
z=interpLagran(orbit(j-5:j+5,1),orbit(j-5:j+5,4),time);
dt=interpLagran(orbit(j-5:j+5,1),orbit(j-5:j+5,5),time);
pos=[x y z]';
dt=dt*1e-6;
end