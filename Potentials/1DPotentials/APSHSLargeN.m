function [VBack_S,VAdd_S]=APSHSLargeN(y,t,optsPhys)


    Vm      = optsPhys.Vm;
    yIh     = optsPhys.yIh;
    Z       = optsPhys.Z;
    
    b=30;
    a=100;
    yI=2;

    if(t==0)
        yt = optsPhys.yt;
    else
        yt = optsPhys.yt2;
    end

    
    h       =  (erf((y+yt)/yIh) - erf((y-yt)/yIh))/2;
    Dh      =  1/sqrt(pi)/yIh*( exp(-(y+yt).^2/yIh^2) - exp(-(y-yt).^2/yIh^2) );

    %--------------------------------------------------------------------------
    VBack  = Vm*y.^2;
    DVBack = 2*Vm*y;
    %--------------------------------------------------------------------------

    VAdd    =  Vm*( h.*(b  -y.^2) - a*Z.*exp(-(y-yt).^2 /yI^2 ) ) ;
    VAdd(abs(y)==inf) = 0;
    DVAdd   =  Vm*( Dh.*(b  -y.^2) + h.*(-2*y) ...
                + a* Z.*exp(- ( y-yt).^2 /yI^2 ).*2.*(y-yt)/yI^2 );
    DVAdd(abs(y)==inf) = 0;
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end
