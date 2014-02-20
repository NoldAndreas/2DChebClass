function [muSC,fnCS,dmuSC,ddmuSC] = MuCarnahanStarling(rho,kBT)
%[muSC,fnCS] = MuCarnahanStarling(rho,kBT)
% muSC  = d(rho*f(rho))/d(rho)
% fnCS  = rho*f(rho)
% dmuSC = d(muSC)/drho
    y      = pi/6*rho;
    muSC   = kBT*(y.*(8-9*y+3*y.^2)./(1-y).^3);
    fnCS   = kBT*rho.*(y.*(4-3*y))./((1-y).^2);
    dmuSC  = kBT*(pi/6)*2*(4-y)./((1-y).^4);
    ddmuSC = kBT*(pi/6)^2*6*(5-y)./((1-y).^5);
    
end
