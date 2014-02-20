function [VBack_S,VAdd_S]=Vext_Cart_TestCahnHillard_1(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


    v2struct(optsPhys);

    VBack        = zeros(size(y1));

    DVBackDy1    = zeros(size(y1));
    DVBackDy2    = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);
                        
    epsilon_w = optsPhys.epsilon_w;
    
    if(~isfield(optsPhys,'epsilon_w_end'))
        epsilon_w_end = 2*epsilon_w;
    else
        epsilon_w_end = optsPhys.epsilon_w_end;
    end
    
    VAdd    = y2*(epsilon_w + (epsilon_w_end - epsilon_w)*(1 -exp(-t^2)));
    VAdd_S  = struct('V',VAdd);

end