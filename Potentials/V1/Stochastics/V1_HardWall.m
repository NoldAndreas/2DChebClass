function [VGeom_S] = V1_HardWall(y1S,y2S,t,optsPhys)

kBT = 1;

if(isfield(optsPhys,'sigma'))
    sigma = optsPhys.sigma;
else
    sigma = 0;
end

%--------------------------------------------------------------------------

V=zeros(size(y1S));
DV1=V;
DV2=V;

% set potential for walls

V1Cutoff = 1/2*sigma;

maskY = (y2S < V1Cutoff);

V(maskY)=NaN;

%--------------------------------------------------------------------------

% smoothed gradient of potential
% V = (sigma/2/r)^2*n - (sigma/2/r)^n + 1/4  if r<2^1/n sigma/2
%     0                                  else       

% V' = (-2*n)*(sigma/2)^(2*n) * (1/r)^(2*n+1) ...
%          + n*(sigma/2)^n * (1/r)^(n+1)           if r<2^1/n sigma/2
% DV = V' dr/dx = x/r V'

% acts perpendicularly to the walls so get
% x/r = 1  for X0, Y0
% x/r = -1 for XL, YL

n=6;

if (sigma~=0)
    sigmaDV = sigma;
else
    sigmaDV = 1;
end

DV1Cutoff = 2^(1/n)*sigmaDV/2;

maskY0 = (y2S < DV1Cutoff);


DV2(maskY0) = DV2(maskY0) + kBT*( - 2*n*(sigmaDV/2)^(2*n) .*y2S(maskY0).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*y2S(maskY0).^(-(n+1))  );
                                                                      
%--------------------------------------------------------------------------
                             
VGeom_S = struct('V',V, ...
                    'dy1',DV1,'dy2',DV2, ...
                    'grad', [DV1;DV2]);                             
                             
end