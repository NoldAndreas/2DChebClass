function [VBack_S,VAdd_S,VGeom_S] = quad_Box(y1S,y2S,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
y1S=y1S(:);

tau = optsPhys.tau;
y10 = optsPhys.y10;
y20 = optsPhys.y20;
B10 = optsPhys.B10;
B20 = optsPhys.B20;
V0Add  = optsPhys.V0Add;

VBack=zeros(size(y1S));
DVBack=zeros(size(y1S));

t = t/tau;

V0Add = V0Add * exp(-t);

VAdd        = - V0Add.*exp(-(y1S-y10).^2/B10).*exp(-(y2S-y20).^2/B20);

DVAddDy1    = - 2*VAdd.*(y1S-y10)/B10;
DVAddDy2    = - 2*VAdd.*(y2S-y20)/B20;
 
DVAddDy1(abs(y1S)==inf)=0;
DVAddDy2(abs(y2S)==inf)=0;

%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DVBack,'dy2',DVBack, ...
                    'grad', [DVBack;DVBack]);

VAdd_S  = struct('V',VAdd, ...
            'dy1',DVAddDy1,'dy2',DVAddDy2,...
            'grad',[DVAddDy1;DVAddDy2]);

otpsPhys.kBT = 1;
        
VGeom_S = V1_Box(y1S,y2S,t,optsPhys);
                
end
