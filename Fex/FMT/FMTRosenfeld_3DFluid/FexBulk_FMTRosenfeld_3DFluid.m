function [mu,fn,dmu,ddmu] = FexBulk_FMTRosenfeld_3DFluid(rho,kBT)
    y   = pi/6*rho;
    
    %here, we compute dOmega/drho    
    mu   = kBT*(-log(1-y)+y.*(14-13*y+5*y.^2)./(2*(1-y).^3));    
    dmu  = kBT*(-(-8+2*y-4*y.^2+y.^3)./(-1+y).^4)*pi/6;
    fn   = kBT*rho.*(-log(1-y) + y*(6-3*y)./(2*(1-y).^2));
    ddmu = kBT*(pi/6)^2*(30+2*y+5*y.^2-y.^3)./((1-y).^5);
end