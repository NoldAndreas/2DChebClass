function [VBack_S,VAdd_S]=EricV1(y,t,optsPhys)


    Vm      = optsPhys.Vm;
    Vq      = optsPhys.Vq;
    alpha   = optsPhys.alpha;
    yt      = optsPhys.y0;

    h       =  (erf((y+yt)/alpha) - erf((y-yt)/alpha))/2;
    Dh      =  1/sqrt(pi)/alpha*( exp(-(y+yt).^2/alpha^2) - exp(-(y-yt).^2/alpha^2) );
    
    %--------------------------------------------------------------------------
    VBack  = Vm.*y.^2;
    DVBack = 2*Vm.*y;
    %--------------------------------------------------------------------------

%     if(t==0)
%         yt = 0;
%     end
    
    VAdd    =  Vq.*h.*(yt - y.^2);
    VAdd(abs(y)==inf) = 0;
    DVAdd   =  Vq.*Dh.*(yt - y.^2) + Vq.*h.*(-2*y);
    DVAdd(abs(y)==inf) = 0;
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end
