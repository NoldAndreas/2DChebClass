function [muSC,fnCS,dmuSC,ddmuSC] = DoublewellPotential(rho,Cn)

% muSC  = d(rho*f(rho))/d(rho)
% fnCS  = rho*f(rho)
% dmuSC = d(muSC)/drho    

    muSC   = -2*rho.*(1-rho.^2)/Cn;
    fnCS   = (1-rho.^2).^2/(2*Cn);        
    dmuSC  = 2*(3*rho.^2-1)/Cn;
    ddmuSC = 12*rho/Cn;
    
end
