function C=DCM_en(L,lon,varargin)
% Compute DCM from navigation to ECEF frame
% input-> lat,lon (deg)
% output-> Cen
if nargin==2
    frame='NED';
else
    frame=varargin{end};
end
switch frame
    case 'NED'
C=[-sind(L)*cosd(lon) -sind(lon) -cosd(L)*cosd(lon);...
    -sind(L)*sind(lon) cosd(lon) -cosd(L)*sind(lon);...
    cosd(L) 0 -sind(L)];
    case 'ENU'
    C=0;
end

end
