function C=rotX(theta)
% Rotate about x-axis following right-handed convention
% input -> theta 
% output -> Rx 
C=[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
end
