 function AddPaths(dirOrg)
   
    global dirData
    global dirDDFT
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
        no = fprintf('Do you want to be asked if data should be recomputed? (press any key)\n');            
        if(getkeywait(2) == -1)
            loadAll = true;
        else
            cprintf('*g','Thanks. You will be asked if data should be recomputed.\n');
            loadAll = false;
        end
    end

    if(isempty(QuickOutput))
        QuickOutput = false;
    end
    
    %recomputeAll = true;
    
    if(recomputeAll)
        cprintf('*m','!!! No precomputed data will be used. recomputeAll = true !!!\n');
    elseif(loadAll)
        cprintf('*m','Data will be used loaded if available. loadAll = true \n');
    end
end