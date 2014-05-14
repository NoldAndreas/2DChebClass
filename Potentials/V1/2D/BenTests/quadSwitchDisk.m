function [VBack_S,VAdd_S]=quadSwitchDisk(r,theta,t,optsPhys)

V0=optsPhys.V0;
y10=optsPhys.y10;
y20=optsPhys.y20;

tau=optsPhys.tau;

%--------------------------------------------------------------------------   

VBack        = zeros(size(r));

DVBackDy1    = zeros(size(r));
DVBackDy2    = zeros(size(r));

VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

%--------------------------------------------------------------------------
        
t = t/tau;

tSwitch = (1 -exp(-t^2));

r(r==inf)=0;

y1S=r.*cos(theta);
y2S=r.*sin(theta);

d = (y1S-y10).^2 + (y2S-y20).^2;

VAdd        = tSwitch*V0.* d;

DVAddPrime = tSwitch*V0;
DdDr     = 2*(y1S-y10).*cos(theta) + 2*(y2S-y20).*sin(theta);
%DdDtheta = -2*(y1S-y10).*r.*sin(theta) + 2*(y2S-y20).*r.*cos(theta);
DdDtheta_r = -2*(y1S-y10).*sin(theta) + 2*(y2S-y20).*cos(theta);

DVAddDr     = DVAddPrime.*DdDr;

DVAddDtheta_r = DVAddPrime.*DdDtheta_r;

VAdd_S  = struct('V',VAdd, ...
            'dy1',DVAddDr,'dy2',DVAddDtheta_r,...
            'grad',[DVAddDr;DVAddDtheta_r]);



end