 function AddPaths()   

    global dirData
    global dirDataOrg
    global dirResults
    global PersonalUserOutput
    global MinimalOutput
    global CodeVersionsDir
    
    global recomputeAll
     
    if(exist('D:\','dir'))
        dirData    = 'D:\2DChebData';    
        dirResults = 'D:\Results\2DCode';        
        dirDDFT    = pwd;
    elseif(exist('/Users/NoldAndreas/','dir'))
        dirData    = '/Users/NoldAndreas/Documents/2DChebData';
        dirResults = '/Users/NoldAndreas/Documents/Results/2DCode/';    
        dirDDFT    = pwd;
    elseif(exist('/home/an2609/','dir'))
        dirData    = '/home/an2609/2DChebData';
        dirDDFT    = pwd;
        dirResults = '/home/an2609/Results';
    elseif(exist('/home/bgoddard/','dir'))
        dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
        dirDDFT    = '/home/bgoddard/work/MATLAB/Fluids/2DChebClass';
        dirResults = '/home/bgoddard/work/MATLAB/Fluids/2DChebData/';
    elseif(exist('/Users/Ben/','dir'))
        dirData    = '/Users/Ben/work/MATLAB/Fluids/2DChebData';
        dirDDFT    = '/Users/Ben/work/MATLAB/Fluids/2DChebClass';
        dirResults = '/Users/Ben/work/MATLAB/Fluids/2DChebData/';
    else
        disp('Unknown computer; using current directory for save data');
        dirData     = pwd;
        dirDDFT     = pwd;
        dirResults  = pwd;
    end
    
    CodeVersionsDir = [dirData filesep 'CodeVersions'];
    
    addpath(genpath(dirDDFT));        
    rmpath(genpath([pwd filesep 'NoClass']));       
    
    PersonalUserOutput = true;
    MinimalOutput      = false;
    
    dirDataOrg = dirData;
    
    recomputeAll = true;
     
end