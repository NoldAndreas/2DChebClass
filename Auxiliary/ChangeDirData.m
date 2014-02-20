function ChangeDirData(newDirData,ORG)

    global dirData
    global dirDataOrg
    
    if( (nargin == 2) && strcmp(ORG,'ORG'))        
        dirDataOrg = newDirData;        
        dirData    = dirDataOrg;

        if(~exist(dirDataOrg,'dir'))            
            disp('Folder not found. Creating new path for dirData..');            
            mkdir(dirDataOrg);
        end
    else
        if(nargin == 0)        
            dirData = dirDataOrg;
        else
            dirData = newDirData;
        end

        if(~exist(dirData,'dir'))            
            disp('Folder not found. Creating new path for dirData..');            
            mkdir(dirData);
        end
    end
    

end