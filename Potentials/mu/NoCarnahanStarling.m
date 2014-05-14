function [muSC,fnCS] = NoCarnahanStarling(rho,kBT)
%     y    = pi/6*rho;
%     muSC = kBT*(y.*(8-9*y+3*y.^2)./(1-y).^3);
%     fnCS = kBT*rho.*(y.*(4-3*y))./((1-y).^2);

muSC=zeros(size(rho));
fnCS=zeros(size(rho));

end
