function [VBack_S,VAdd_S]=quadQuadBump(y1S,y2S,t,optsPhys)

V0x=optsPhys.V0x;
V0y=optsPhys.V0y;
V0Add = optsPhys.V0Add;
A0 = optsPhys.A0;
B0 = optsPhys.B0;
tau=optsPhys.tau;

%--------------------------------------------------------------------------   

VBack        = V0x.*y1S.^2 + V0y.*y2S.^2;

DVBackDy1    = 2*V0x.*y1S;
DVBackDy2    = 2*V0y.*y2S;

DVBackDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
DVBackDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;

VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

%--------------------------------------------------------------------------
        
t = exp(-t/tau);

V0Add = t*V0Add;

VAdd        = V0Add.*exp(-(y1S-A0).^2/B0);

DVAddDy1    = -2*VAdd.*(y1S-A0)/B0;
DVAddDy2    = zeros(size(y2S));
 
DVAddDy1(abs(y1S)==inf)=0;
DVAddDy2(abs(y2S)==inf)=0;

VAdd_S  = struct('V',VAdd, ...
            'dy1',DVAddDy1,'dy2',DVAddDy2,...
            'grad',[DVAddDy1;DVAddDy2]);



end