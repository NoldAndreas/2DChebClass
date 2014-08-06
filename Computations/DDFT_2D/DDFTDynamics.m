function EX = DDFTDynamics(optsPhys,optsNum,optsPlot)

    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();

    if(EX.doHIWall)
        EX.ComputeDynamicsWallHI();
    elseif(isfield(optsPhys,'Inertial') && optsPhys.Inertial)
        EX.ComputeDynamicsInertia();
    else
        EX.ComputeDynamics();
    end
   
    if(isfield(optsPlot,'doDDFTPlots') && optsPlot.doDDFTPlots)
        EX.PlotDynamics();
    end
    
    if(isfield(optsNum,'doPlots') && optsNum.doPlots)
        EX.PlotDynamics();
    end
end