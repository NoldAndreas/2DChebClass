function [mu_s,fnCS,dmuSC,ddmuSC] = FexBulk_CarnahanStarling(rho_s,kBT)
    [mu_s,fnCS,dmuSC,ddmuSC]  = CarnahanStarling(rho_s,kBT);
end