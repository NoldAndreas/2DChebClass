 function AddPaths(dirOrg)
   
    global dirData
    global dirDataOrg    
    
    global PersonalUserOutput
    global QuickOutput
    
    global recomputeAll
    global loadAll
     
    if(exist('D:\','dir'))
        dirData    = 'D:\2DChebData';    
        dirDDFT    = pwd;
    elseif(exist('/Users/NoldAndreas/','dir'))
        dirData    = '/Users/NoldAndreas/Documents/2DChebData';
        dirDDFT    = pwd;
    elseif(exist('/home/an2609/','dir'))
        dirData    = '/home/an2609/2DChebData';
        dirDDFT    = pwd;        
    elseif(exist('/home/bgoddard/','dir'))
        dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
        dirDDFT    = '/home/bgoddard/work/MATLAB/Fluids/2DChebClass';        
    elseif(exist('/Users/Ben/','dir'))
        dirData    = '/Users/Ben/work/MATLAB/Fluids/2DChebData';
        dirDDFT    = '/Users/Ben/work/MATLAB/Fluids/2DChebClass';        
    else
        disp('Unknown computer; using current directory to save data');
        dirData     = pwd;
        dirDDFT     = pwd;        
    end
    
    addpath(genpath(dirDDFT));        
    rmpath(genpath([pwd filesep 'NoClass']));       
    
    PersonalUserOutput = true;    
    
    dirDataOrg = dirData;
    if(nargin >= 1)
        ChangeDirData([dirDataOrg filesep dirOrg],'ORG');
    end
    
    if(isempty(recomputeAll))
        recomputeAll = false;
    end
    
	if(isempty(loadAll))
        loadAll = false;
    end
    
    if(isempty(QuickOutput))
        QuickOutput = false;
    end        
    
    if(recomputeAll)
        cprintf('*m','!!! No precomputed data will be used. recomputeAll = true !!!\n');
    end
end