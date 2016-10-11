function [VBack_S,VAdd_S]=gravity3D(y1S,y2S,y3S,t,optsPhys)

    g = optsPhys.g;

    %--------------------------------------------------------------------------
    VBack           = zeros(size(y1S));
    
    DVBackDy1          = zeros(size(y1S));        
    DVBackDy2          = zeros(size(y1S));        
    DVBackDy3          = zeros(size(y1S));        
    %--------------------------------------------------------------------------
    
    VAdd          = g*y3S;
    
    DVAddDy1     = zeros(size(y1S));
    DVAddDy2     = zeros(size(y1S));
    DVAddDy3     = g*ones(size(y1S));
        
    %--------------------------------------------------------------------------
    
    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,'dy3',DVBackDy3, ...
                'grad',[DVBackDy1;DVBackDy2;DVBackDy3]);

    %--------------------------------------------------------------------------
        

    VAdd_S  = struct('V',VAdd, ...
                'dy1',DVAddDy1,'dy2',DVAddDy2,'dy3',DVAddDy3,...
                'grad',[DVAddDy1;DVAddDy2;DVAddDy3]);
        
end