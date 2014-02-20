function [mu_s,fnCS,dmuSC,ddmuSC] = FexBulk_CarnahanStarling(rho_s,kBT)
    [mu_s,fnCS,dmuSC,ddmuSC]  = MuCarnahanStarling(rho_s,kBT);
end