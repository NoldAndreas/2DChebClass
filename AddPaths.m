 function AddPaths(dirOrg)
   
    global dirData
    global dirDDFT
    global dirDataOrg    
    
    global PersonalUserOutput
    global QuickOutput
    
    global recomputeAll
    global loadAll         
    
    switch GetMacAddress()
        case '24-BE-05-10-A1-52'  %Andreas' Windows Work PC
            dirData    = 'D:\2DChebData';    
            dirDDFT    = pwd;
        case '00:88:65:35:a1:92'
            dirData    = '/Volumes/BACKUP_IC/2DChebData';
            %dirData    = '/Users/NoldAndreas/Documents/2DChebData';
            dirDDFT    = pwd;
        case '1C:C1:DE:52:03:FF' % Ben office machine
            dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
            dirDDFT    = pwd;
        case '10:40:f3:8a:30:f4' % Ben MacBook Air
            dirData    = '/Users/Ben/work/MATLAB/Fluids/2DChebData';
            dirDDFT    = pwd;
        case 'C8:1F:66:ED:06:8F' % compute64c
            switch getenv('USER')
                case 'bgoddard'
                    dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
                    dirDDFT    = pwd;
                otherwise
                    disp('Unknown computer; using current directory to save data');
                    dirData     = pwd;
                    dirDDFT     = pwd;        
            end
        case 'A0:36:9F:60:6C:C4' % compute64d
            switch getenv('USER')
                case 'bgoddard'
                    dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
                    dirDDFT    = pwd;
                otherwise
                    disp('Unknown computer; using current directory to save data');
                    dirData     = pwd;
                    dirDDFT     = pwd;        
            end
        otherwise
            disp('Unknown computer; using current directory to save data');
            dirData     = pwd;
            dirDDFT     = pwd;        
    end

%    elseif(exist('/home/an2609/','dir'))
%        dirData    = '/home/an2609/2DChebData';
%        dirDDFT    = pwd;
% %     elseif(exist('/home/bgoddard/','dir'))
% %         dirData    = '/home/bgoddard/work/MATLAB/Fluids/2DChebData';
% %         dirDDFT    = '/home/bgoddard/work/MATLAB/Fluids/2DChebClass';        
% %     elseif(exist('/Users/Ben/','dir'))
% %         dirData    = '/Users/Ben/work/MATLAB/Fluids/2DChebData';
% %         dirDDFT    = '/Users/Ben/work/MATLAB/Fluids/2DChebClass';                    
% %     end
    
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
    
    function mac_add = GetMacAddress()
        switch computer('arch')
             case {'maci','maci64'}
                 [~,a]=system('ifconfig');
                 c=strfind(a,'en0');if ~isempty(c),a=a(c:end);end
                 c=strfind(a,'en1');if ~isempty(c),a=a(1:c-1);end
                 % find the mac address
                 b=strfind(a,'ether');
                 mac_add=a(1,b(1)+6:b(1)+22);
             case {'win32','win64'}
                 [~,a]=system('getmac');
                 b=strfind(a,'=');
                 %mac_add=a(b(end):b(end)+19);
                 mac_add=a(b(end)+2:b(end)+18);
             case {'glnx86','glnxa64'}
                 [~,a]=system('ifconfig');b=strfind(a,'Ether');
                 mac_add=a(1,b(1)+17:b(1)+33);
             otherwise,mac_add=[];
        end

        mac_add = strtrim(mac_add);
    end
end