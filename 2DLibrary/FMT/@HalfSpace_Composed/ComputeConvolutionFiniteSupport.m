function [AAD] = ComputeConvolutionFiniteSupport(this,area,weights,pts)
%GetAverageDensities(this,area,weights,pts)            
    %AAD - average the average densities to compute free energy
    fprintf('Computing interpolation for matrices for averaging of averaged densities..\n');    
    AAD  = Conv_LinearGridXY(this,pts,area,weights);        
end  