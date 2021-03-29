function [n,e,d] = nedFromEcef(x,y,z,lat0,lon0)
if size(x,1)==1
    x=x';
    y=y';
    z=z';
end
p=[x-x(1) y-y(1) z-z(1)];
Cen=DCM_en(lat0,lon0);
p=Cen*p';
n=p(1,:);
e=p(2,:);
d=p(3,:);
end

