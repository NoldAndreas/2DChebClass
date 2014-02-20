function [mu_s,h,fnCS] = Fex_CarnahanStarling(rho_s,IntMatrFex,kBT,R)
    [mu_s,fnCS,dmuSC,ddmuSC]  = MuCarnahanStarling(rho_s,kBT);
    h = [];
end