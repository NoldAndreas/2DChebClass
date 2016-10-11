function [VBack_S,VAdd_S]=EricV1(y,t,optsPhys)


    Vm      = optsPhys.Vm;
    Vq      = optsPhys.Vq;
    alpha   = optsPhys.alpha;
    
    if(t==0)
        yt = optsPhys.y0;
    else
        yt = optsPhys.y1;
    end    

    h       =  (erf((y+yt)/alpha) - erf((y-yt)/alpha))/2;
    Dh      =  1/sqrt(pi)/alpha*( exp(-(y+yt).^2/alpha^2) - exp(-(y-yt).^2/alpha^2) );
    
    %--------------------------------------------------------------------------
    VBack  = Vm.*y.^2;
    DVBack = 2*Vm.*y;
    %--------------------------------------------------------------------------
    
    VAdd    =  Vq.*h.*(yt^2 - y.^2);
    VAdd(abs(y)==inf) = 0;
    DVAdd   =  Vq.*Dh.*(yt^2 - y.^2) + Vq.*h.*(-2*y);
    DVAdd(abs(y)==inf) = 0;
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end
