function [VBack_S,VAdd_S]=Vext_Cart_Capillary_2(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0              = optsPhys.V0;
    y10             = optsPhys.y10;
    y20             = optsPhys.y20;
    t               = t/optsPhys.tau;
    
    yw_left   = min(y1); %yw1
    yw_right  = max(y1); %yw1
    yw_bottom = min(y2); %yw2
    yw_top    = max(y2); %yw3            
        
    VBack    = V0*((y1-y10).^2); %+ (y2-y20).^2);
    
    DVBackDy1 = 2*V0*(y1-y10);% + da1;
    DVBackDy2 = zeros(size(y1));%2*V0*(y2-y20);%+ da2 - da3;

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            

    epsilon_w = optsPhys.epsilon_w + (1 - exp(-t.^2))*...
                        (optsPhys.epsilon_w_end - optsPhys.epsilon_w);                    
    
    d      = yw_top - yw_bottom;
    a1     = Slit(y2-yw_bottom,yw_right-y1,d,epsilon_w(1));
    a2     = AttractiveWall(yw_top-y2,epsilon_w(2),1);
    a3     = Slit(y2-yw_bottom,y1 - yw_left,d,epsilon_w(3));    
    a4     = AttractiveWall(y2-yw_bottom,epsilon_w(4),1);    
    
    b      = -(1-exp(-t^2))*V0*((y1-y10).^2);
    
    VAdd    = a1 + a2 + a3 + a4 + b;
    VAdd_S  = struct('V',VAdd);

end