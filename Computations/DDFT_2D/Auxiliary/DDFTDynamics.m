function EX = DDFTDynamics(optsPhys,optsNum,optsPlot)

    AddPaths();
    EX   = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();
 
    if(EX.doHIWall)
        EX.ComputeDynamicsWallHI();
    elseif( (isfield(optsPhys,'Inertial') && optsPhys.Inertial) || ...
             (isfield(optsNum,'Inertial') && optsNum.Inertial) )
        EX.ComputeDynamicsInertia();
    else
        EX.ComputeDynamics();
    end
   
    if( (nargin < 3) || ...
        (isfield(optsPlot,'doDDFTPlots') && optsPlot.doDDFTPlots) || ...
        (isfield(optsNum,'doPlots') && optsNum.doPlots) )
        EX.PlotDynamics();
    end
    
end