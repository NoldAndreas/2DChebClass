function [y1R,y2R,J] = RotateAroundOrigin(y1,y2,th)
    y1R = cos(th)*y1 - sin(th)*y2;
	y2R = sin(th)*y1 + cos(th)*y2;
    
    n        = length(y1);
	J        = zeros(n,2,2);
    J(:,1,1) = cos(th);
    J(:,1,2) = -sin(th);
    J(:,2,1) = sin(th);
    J(:,2,2) = cos(th);
    
end
