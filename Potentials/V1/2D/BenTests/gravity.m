function [VBack_S,VAdd_S] = gravity(y1S,y2S,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
y1S=y1S(:);

g     = optsPhys.g;
if(isfield(optsPhys,'theta'))
    theta = optsPhys.theta;
else
    theta = 0;
end

theta = theta*pi/180;

VBack  = zeros(size(y1S));
DVBack = zeros(size(y1S));


g = g * (1- exp(-t));

h = sin(theta).*y1S + cos(theta).*y2S;

VAdd        = g*h;

DVAddDy1    = g*sin(theta)*ones(size(y1S));
DVAddDy2    = g*cos(theta)*ones(size(y2S));
 
%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DVBack,'dy2',DVBack, ...
                    'grad', [DVBack;DVBack]);

VAdd_S  = struct('V',VAdd, ...
            'dy1',DVAddDy1,'dy2',DVAddDy2,...
            'grad',[DVAddDy1;DVAddDy2]);
                
end
