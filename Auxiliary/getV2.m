function [z,dzdr_r,alpha] = getV2(r,V2)

global r_cutoff
if(isfield(V2,'r_cutoff'))    
    r_cutoff = V2.r_cutoff;
end

% get name of function
if(ischar(V2.V2DV2))        
    getV2DV2 =str2func(V2.V2DV2);
else
    getV2DV2 = V2.V2DV2;
end

% retrieve rescaled values
if(nargin(getV2DV2) == 1)
    [z,dzdr_r,alpha] = getV2DV2(r);

    z      = V2.epsilon*z;
    dzdr_r = V2.epsilon*dzdr_r;
    alpha  = V2.epsilon*alpha;
else
    [z,dzdr_r,alpha] = getV2DV2(r,V2);
end
    
