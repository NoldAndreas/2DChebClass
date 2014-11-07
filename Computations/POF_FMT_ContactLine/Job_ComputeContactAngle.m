function Job_ComputeContactAngle(opts)
    
    config = GetStandardConfig(opts);
    
    %***********************************************************        
    %filename   = [dirData,filesep,'Job40_DisjoiningPressure_',getTimeStr(),'.txt'];       
    %Struct2File(filename,config,['Computed at ',datestr(now)]);

  %  ConvergenceSurfaceTensions(config);
    close all;
    
    CLT = ContactLineHS(config);     
    CLT.Preprocess();
    CLT.ComputeEquilibrium();          	        
    CLT.PostProcess();
    
    CLT.PlotDensityResult();
    %CLT.PlotContourResults();
	%CLT.PlotDisjoiningPressures();

    %     CLT.ComputeAdsorptionIsotherm('load'); %load \2DChebData\POF_FMT_ContactLine\deg90\IterativeContinuationPostProcess\2014_1_20_18_46
%     CLT.PostProcess_2DDisjoiningPressure();
%     
% 	%f1 = figure('Color','white','Position',[0 0 1000 600]);
%     [f1,f2] = CLT.Post_HFrom2DDisjoiningPressure();
%
%     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces'],f1);
%     saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces.fig']);
%
%     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure'],f2);
%     saveas(f2,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure.fig']);
%     
%     inset2(f1,f2,0.4,[0.26,0.55]);
%     %inset2(f1,f2,0.35,[0.22,0.55]);
%     close(f2);      
%     
%     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure'],f1);
%     saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure.fig']);
%         
%     CLT.FittingAdsorptionIsotherm([10 14],1);
%     CLT.SumRule_DisjoiningPotential();
%     %***************************************  

end