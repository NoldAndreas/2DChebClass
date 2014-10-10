function TestAll(dirTest)
    global recomputeAll QuickOutput dirData dirDataOrg

    recomputeAll = true;
    QuickOutput  = true;
    comments = 'All computation times given in seconds.';
        
    
    dirResOld  = dirData;
    dirData    = [dirDataOrg,dirTest,'\'];
    
    if(~exist(dirData,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(dirData);
   end
    

    dirDDFT_2D  = [pwd,'/Computations/',dirTest];
    MFiles      = dir(fullfile(dirDDFT_2D,'*.m'));

    Parameters = LoadIndexFiles(dirDDFT_2D);
    if(isempty(Parameters))
        Parameters = struct();
    else
        Parameters = Parameters{1};
    end
    
    for i = 1:length(MFiles)
        strf = MFiles(i).name;
        strf = strf(1:end-2);
        if(~isfield(Parameters,strf))
           Parameters.(strf) = num2str(-1);
        end
    end
    
    no = fprintf('Do you want to rerun All files? (press any key), or wait for 2 seconds.');        
    if(getkeywait(2) == -1)
        for ih = 1:no
            fprintf('\b');
        end                        
        rerun = false;
    else
        for ih = 1:no
            fprintf('\b');
        end
        no = fprintf('Thanks. All files will be rerun.\n');
        rerun = true;
    end                

    %Run through all files
    for i = 1:length(MFiles)
        strf = MFiles(i).name;
        strf = strf(1:end-2);

        if(rerun ||...
           ~isfield(Parameters,strf) || ...
           (str2double(Parameters.(strf)) < 0))
           
            tStart = tic;
            cprintf('*Blue',['** Running ',strf,' **\n']);            
            f = str2func(strf);
            
            try                
                f();            
                Parameters.(strf) = num2str(round(toc(tStart)));
            catch err
                cprintf('*Red',['Error in ',strf,'\n']);
                Parameters.(strf) = num2str(-2);
            end
            close all;
        else
            cprintf('*Blue',['** Skipping ',strf,' **\n']);
        end
        
        WriteReport();
    end
    
    dirData = dirResOld;

    function WriteReport()
        tc = 0;
        for j = 1:length(MFiles)
            strfT = MFiles(j).name;
            strfT = strfT(1:end-2);
            if(str2num(Parameters.(strfT)) > 0)
               tc = tc + str2num(Parameters.(strfT));
            end
        end

        Struct2File([dirDDFT_2D,'/TestAll_report.txt'],Parameters,...
                        ['Computed at: ',datestr(now),...
                        ' Cumulative computation time: ',sec2hms(tc), ' (hrs:min:sec) ',comments]);

    end
                    
    recomputeAll = false;
end
