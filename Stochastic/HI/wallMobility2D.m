function [D,sumDivD] = wallMobility2D(r,optsPhys)

% See von Hansen, Hinczewski, Netz, JCP, 134, 235102 (2011)
% only works for one species

D0=optsPhys.D0;
sigmaH=optsPhys.sigmaH;

sigma = sigmaH(1);

%--------------------------------------------------------------------------
% Extract x and z components 
%--------------------------------------------------------------------------

x = r(1:2:end);
z = r(2:2:end);

% Self mobility terms are included automatically (except D0 term from
% identity part of HI)
zInv = 1./z;
selfTermPar = - 9/32*sigma*zInv + 1/64*(sigma*zInv).^3;
selfTermPer = - 9/16*sigma*zInv + 1/16*(sigma*zInv).^3;

Dm = [selfTermPar'; selfTermPer'];
Dm = Dm(:);
Dm = diag(Dm);

% calculate Rotne-Prager form
D=diag(D0)*(eye(length(r)) + Dm);

divPar = zeros(size(x));
divPer = 9*sigma/32 * zInv.^2 - 3*sigma^3/16 * zInv.^4;

sumDivD = [divPar';divPer'];
sumDivD = sumDivD(:);

