function C=rotY(theta)
% Rotate about y-axis following right-handed convention
% input -> theta 
% output -> Ry
C=[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
end
