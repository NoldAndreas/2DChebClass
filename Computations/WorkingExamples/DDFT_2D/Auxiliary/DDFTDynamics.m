function [EX,res] = DDFTDynamics(optsPhys,optsNum,optsPlot)

    AddPaths();
    EX   = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    

    EX.ComputeEquilibrium([],optsNum.eqNum); 
    
    EX.ComputeDynamics();
   
    if( (nargin < 3) || ...
        (isfield(optsPlot,'doDDFTPlots') && optsPlot.doDDFTPlots) || ...
        (isfield(optsNum,'doPlots') && optsNum.doPlots) )
        res.fig_handles = EX.PlotDynamics();
    end
    
end