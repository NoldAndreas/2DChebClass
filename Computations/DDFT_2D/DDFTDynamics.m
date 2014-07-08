function EX = DDFTDynamics(optsPhys,optsNum,optsPlot)

    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();
    
    if(EX.doHIWall)
        EX.ComputeDynamicsWallHI();
    else
        EX.ComputeDynamics();
    end
    
    if(optsPlot.doDDFTPlots)
        EX.PlotDynamics();
    end
end