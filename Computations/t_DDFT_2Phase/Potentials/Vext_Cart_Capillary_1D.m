function [VBack_S,VAdd_S]=Vext_Cart_Capillary_1D(y,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    %V0              = optsPhys.V0;
    %y10             = optsPhys.y10;
    %y20             = optsPhys.y20; 
    t               = t/optsPhys.tau;
    
    yw2 = min(y);
    yw3 = max(y);        
        
    VBack    = zeros(size(y));%V0*((y1-y10).^2); %+ (y2-y20).^2);
    
    DVBackDy = zeros(size(y));%2*V0*(y1-y10);% + da1;    

    VBack_S = struct('V',VBack,'dy',DVBackDy);            

    epsilon_w2 = optsPhys.epsilon_w2 + (1 - exp(-t.^2))*...
                        (optsPhys.epsilon_w2_end - optsPhys.epsilon_w2);
    epsilon_w3 = optsPhys.epsilon_w3 + (1 - exp(-t.^2))*...
                        (optsPhys.epsilon_w3_end - optsPhys.epsilon_w3);                    
                
    [a2,da2] = AttractiveWall(y-yw2,epsilon_w2,1);
    [a3,da3] = AttractiveWall(yw3-y,epsilon_w3,1);
                           
    VAdd    = a2 + a3;
    VAdd_S  = struct('V',VAdd);

end