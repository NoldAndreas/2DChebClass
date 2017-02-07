function Compute(this)
    
    close all;   
    
%    diaryFile = [dirData filesep 'LogFile.txt'];    diary(diaryFile); diary on;

    Preprocess(this);        	
    if(this.optsNum.maxComp_y2 < 0)
        return;
    end
     
    ComputeEquilibrium(this);           
    
    PostProcess(this,[-1,5]);
    PlotEquilibriumResults(this);
    PlotDisjoiningPressureAnalysis(this);    
    
    if(~isfield(this.optsNum,'plotTimes_T'))%diary                
        return;
    end    
    

end