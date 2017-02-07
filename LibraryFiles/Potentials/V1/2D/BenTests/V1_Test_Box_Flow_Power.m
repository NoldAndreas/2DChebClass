function [VBack_S,VAdd_S,VGeom_S]=V1_Test_Box_Flow_Power(y1S,y2S,t,optsPhys)

V0=optsPhys.V0;
y10=optsPhys.y10;
y20=optsPhys.y20;

L1 = optsPhys.L1;
L2 = optsPhys.L2;

epsilon = optsPhys.epsilon0;
alpha   = optsPhys.alphaWall;

pow = optsPhys.pow;

tau=optsPhys.tau;

%--------------------------------------------------------------------------   

VWall11      = epsilon.*exp(-(y1S.^2)./(alpha.^2));
VWall12      = epsilon.*exp(-((y1S-L1).^2)./(alpha.^2));

VWall21      = epsilon.*exp(-(y2S.^2)./(alpha.^2));
VWall22      = epsilon.*exp(-((y2S-L2).^2)./(alpha.^2));

DVWall11     = -2*y1S./(alpha.^2).*VWall11;
DVWall12     = -2*(y1S-L1)./(alpha.^2).*VWall12;

DVWall21     = -2*y2S./(alpha.^2).*VWall21;
DVWall22     = -2*(y2S-L2)./(alpha.^2).*VWall22;


% VWall11 = epsilon.*(y1S./alpha).^(-1);
% VWall12 = epsilon.*((L1-y1S)./alpha).^(-1);
% 
% VWall21 = epsilon.*(y2S./alpha).^(-1);
% VWall22 = epsilon.*((L2-y2S)./alpha).^(-1);
% 
% 
% DVWall11 = epsilon.*(y1S).^(-2).*alpha;
% DVWall12 = -epsilon.*(L1-y1S).^(-2).*alpha;
% 
% DVWall21 = epsilon.*(y2S).^(-2).*alpha;
% DVWall22 = -epsilon.*(L2-y2S).^(-2).*alpha;


VBack        =  VWall11 + VWall12 + VWall21 + VWall22;

DVBackDy1    = DVWall11 + DVWall12;
DVBackDy2    = DVWall21 + DVWall22;

VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

%--------------------------------------------------------------------------
        
t = t/tau;

tSwitch = exp(-t^2);

VAdd        = tSwitch*V0.*((y1S-y10).^pow);

VAdd(abs(y1S)==inf | abs(y2S)==inf) = 0;

DVAddDy1    = tSwitch*pow*V0.*(y1S-y10).^(pow-1);
DVAddDy2    = zeros(size(y2S));

DVAddDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
DVAddDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;

VAdd_S  = struct('V',VAdd, ...
            'dy1',DVAddDy1,'dy2',DVAddDy2,...
            'grad',[DVAddDy1;DVAddDy2]);

%--------------------------------------------------------------------------

VGeom_S = V1_Box(y1S,y2S,t,optsPhys);

end