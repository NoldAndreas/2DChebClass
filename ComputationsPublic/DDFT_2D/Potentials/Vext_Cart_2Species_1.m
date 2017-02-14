function [VBack_S,VAdd_S]=Vext_Cart_2Species_1(y1S,y2S,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


    v2struct(optsPhys);

    VBack        = V0*((y1S-y10).^2 + (y2S-y20).^2);

    DVBackDy1    = 2*V0*(y1S-y10);
    DVBackDy2    = 2*V0*(y2S-y20);

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);
            
    t = t/tau;
    %1st Species
%     y1        = y1S(:,1);
%     wall1     = pi*(atan(y1)-pi/2+y1./(1+y1.^2));
%     VAdd(:,1) = wall1*(epsilon_w1 + (epsilon_w1_end - epsilon_w1)*(1 -exp(-t^2)));
     
%     y2 = y2S(:,1);
%     wall2      = pi*(atan(y2)-pi/2+y2./(1+y2.^2));
%     VAdd(:,2)  = wall2*(epsilon_w2 + (epsilon_w2_end - epsilon_w2)*(1 -exp(-t^2)));
%     VAdd(:,2)  = VAdd(:,2) + wall1*(epsilon_w2_end + (epsilon_w2 - epsilon_w2_end)*(1 -exp(-t^2)));

    wall1     = pi*(atan(y1S)-pi/2+y1S./(1+y1S.^2));
    wall2     = pi*(atan(y2S)-pi/2+y2S./(1+y2S.^2));
    
    
    VAdd      = wall1.*(epsilon_w1 + (epsilon_w1_end - epsilon_w1)*(1 -exp(-t^2)))+ ...
               +wall2.*(epsilon_w2 + (epsilon_w2_end - epsilon_w2)*(1 -exp(-t^2)));    

    VAdd_S  = struct('V',VAdd);

end