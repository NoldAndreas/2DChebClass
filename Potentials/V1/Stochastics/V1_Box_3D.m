function [VGeom_S] = V1_Box_3D(y1S,y2S,y3S,t,optsPhys)

% length of box in each direction
L1    = optsPhys.L1;
L2    = optsPhys.L2;
L3    = optsPhys.L3;
%kBT   = optsPhys.kBT;

kBT = 1;

if(isfield(optsPhys,'sigma'))
    sigma = optsPhys.sigma;
    % sigma is nParticles x nParticles
    sigma = sigma(1);
else
    sigma = 0;
end

%--------------------------------------------------------------------------

V=zeros(size(y1S));
DV1=V;
DV2=V;
DV3=V;

% set potential for walls

V1Cutoff = 1/2*sigma;

maskX0 = (y1S < V1Cutoff);
maskY0 = (y2S < V1Cutoff);
maskZ0 = (y3S < V1Cutoff);

maskXL = (y1S > L1-V1Cutoff);
maskYL = (y2S > L2-V1Cutoff);
maskZL = (y3S > L3-V1Cutoff);

maskAny = maskX0 | maskY0 | maskZ0 | maskXL | maskYL | maskZL;

V(maskAny)=NaN;

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

maskX0 = (y1S < DV1Cutoff);
maskY0 = (y2S < DV1Cutoff);
maskZ0 = (y3S < DV1Cutoff);

maskXL = (y1S > L1-DV1Cutoff);
maskYL = (y2S > L2-DV1Cutoff);
maskZL = (y3S > L3-DV1Cutoff);

DV1(maskX0) = DV1(maskX0) + kBT*( - 2*n*(sigmaDV/2)^(2*n) .*y1S(maskX0).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*y1S(maskX0).^(-(n+1))  );

DV1(maskXL) = DV1(maskXL) - kBT*( - 2*n*(sigmaDV/2)^(2*n) .*(L1-y1S(maskXL)).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*(L1-y1S(maskXL)).^(-(n+1))  );
                                                                                       
DV2(maskY0) = DV2(maskY0) + kBT*( - 2*n*(sigmaDV/2)^(2*n) .*y2S(maskY0).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*y2S(maskY0).^(-(n+1))  );

DV2(maskYL) = DV2(maskYL) - kBT*( - 2*n*(sigmaDV/2)^(2*n) .*(L2-y2S(maskYL)).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*(L2-y2S(maskYL)).^(-(n+1))  );
                             
DV3(maskZ0) = DV2(maskZ0) + kBT*( - 2*n*(sigmaDV/2)^(2*n) .*y2S(maskZ0).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*y2S(maskZ0).^(-(n+1))  );

DV3(maskZL) = DV2(maskZL) - kBT*( - 2*n*(sigmaDV/2)^(2*n) .*(L3-y2S(maskZL)).^(-(2*n+1)) ...
                                 + n*(sigmaDV/2)^n      .*(L3-y2S(maskZL)).^(-(n+1))  );

                             
                             
%--------------------------------------------------------------------------
                             
VGeom_S = struct('V',V, ...
                    'dy1',DV1,'dy2',DV2, 'dy3',DV3, ...
                    'grad', [DV1;DV2;DV3]);                             
                             
end