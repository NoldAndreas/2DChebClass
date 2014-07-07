function EX = DDFTDynamics(optsPhys,optsNum,optsPlot)

    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();
    EX.ComputeDynamics();
    
    if(optsPlot.doDDFTPlots)
        EX.PlotDynamics();
    end
end