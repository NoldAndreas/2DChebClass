function U_S=Flow_Test_Box(y1S,y2S,t,optsPhys)

L2 = optsPhys.BoxHeight;
U0 = optsPhys.U0;

%--------------------------------------------------------------------------   

U1        = U0*y2S.*(L2-y2S);
U2        = zeros(size(y1S));

DU1Dy1    = zeros(size(y1S));
DU2Dy2    = zeros(size(y1S));

U_S = struct('U',[U1;U2],'U1',U1,'U2',U2,'div',DU1Dy1+DU2Dy2);

%--------------------------------------------------------------------------
        
% t = t/tau;
% 
% tSwitch = (1 -exp(-t^2));
% 
% VAdd        = tSwitch*V0.*((y1S-y10).^2 + (y2S-y20).^2);
% 
% VAdd(abs(y1S)==inf | abs(y2S)==inf) = 0;
% 
% DVAddDy1    = 2*V0.*(y1S-y10);
% DVAddDy2    = 2*V0.*(y2S-y20);
% 
% DVAddDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
% DVAddDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;
% 
% VAdd_S  = struct('V',VAdd, ...
%             'dy1',DVAddDy1,'dy2',DVAddDy2,...
%             'grad',[DVAddDy1;DVAddDy2]);

%VGeom_S = V1_Box(y1S,y2S,t,optsPhys);

end