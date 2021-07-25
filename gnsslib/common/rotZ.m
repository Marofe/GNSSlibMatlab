function C=rotZ(theta)
% Rotate about z-axis following right-handed convention
% input -> theta 
% output -> Rz
C=[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
end
