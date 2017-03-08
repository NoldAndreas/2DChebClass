function TestAll(dirTest,recomputeAll,rerun)
    global recomputeAll dirData dirDataOrg
  %  dbclear all;
    
    if((nargin < 2) || isempty(recomputeAll))
        no = fprintf('Do you want to recompute all matrices? (press any key), or wait for 3 seconds.');
        if(getkeywait(3) == -1)
            for ih = 1:no
                fprintf('\b');
            end                        
            recomputeAll = false;
        else
            for ih = 1:no
                fprintf('\b');
            end
            no = fprintf('Thanks. All matrices will be recomputed.\n');
            recomputeAll = true;
        end    
    end
            
    comments     = 'TEST_ALL -- computation time in secs (>0); ignore file (=0); recompute (< 0); error in last run (-2)';
        
    
    dirResOld  = dirData;
    dirData    = [dirDataOrg filesep dirTest filesep];    
    if(~exist(dirData,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(dirData);
    end           

    dirDDFT_2D  = [pwd,filesep,dirTest];
    MFiles      = dir(fullfile(dirDDFT_2D,'*.m'));
    
    
    dirDDFT_2D_LatexReport  = [dirData,'LatexReport'];    
    if(~exist(dirDDFT_2D_LatexReport,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(dirDDFT_2D_LatexReport);
    end           
    filenameLatexDoc = [dirDDFT_2D_LatexReport,'/LatexReport.tex'];
    InitLatexDoc();
    
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
    
    if((nargin < 3) || isempty(rerun))
        no = fprintf('Do you want to rerun All files? (press any key), or wait for 3 seconds.');        
        if(getkeywait(3) == -1)
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
    end

    %Run through all files
    for i = 1:length(MFiles)
        strf = MFiles(i).name;
        strf = strf(1:end-2);
        
        AddToLatexDoc([strf,'.m']);

        if(rerun ||...
           ~isfield(Parameters,strf) || ...
           (str2double(Parameters.(strf)) < 0))
           
            tStart = tic;
            cprintf('*Blue',['** Running ',strf,' **\n']);            
            f = str2func(strf);
            
            try     
                if(nargout(f)>=2)
                    [~,res] = f();
                    
                    if(isstruct(res) && isfield(res,'fig_handles') && isnumeric(res.fig_handles{1}))
                        for k = 1:length(res.fig_handles)                        
                            fh = res.fig_handles{k};                        
                            set(0, 'currentfigure', fh);

                            ax = get(fh,'children');
                            xlim = get(ax,'xlim');
                            ylim = get(ax,'ylim');
                            if(isnumeric(ylim)) %in case of subplots, this is a list
                                r = (ylim(2)-ylim(1))/(xlim(2)-xlim(1));
                                pbaspect(ax,[(xlim(2)-xlim(1)) (ylim(2)-ylim(1)) 1]);
                                
                                hand = get(ax,'xlabel'); set(hand,'str','');                        
                                hand = get(ax,'ylabel'); set(hand,'str','');
                                %hand = get(gca,'xlabel'); set(hand,'fontsize',35);                        
                                %hand = get(gca,'ylabel'); set(hand,'fontsize',35);
                            else
                                r = 1;
                            end
                            set(fh,'Position',[0 0 600 600*r],'color','white');                                                    
                            set(ax,'fontsize',25);                           

                            str_fig = [strf,num2str(k)];
                            set(0, 'currentfigure', fh);
                            SaveFigure([dirDDFT_2D_LatexReport filesep str_fig]);                        
                        end
                    end
                else
                    f();
                end
                Parameters.(strf) = num2str(round(toc(tStart)));
            catch err
                cprintf('*Red',['Error in ',strf,'\n']);         
                Parameters.(strf) = num2str(-2);
            end
            close all;
        else
            cprintf('*Blue',['** Skipping ',strf,' **\n']);
        end
                
        AddToLatexFig([strf,'.eps']);                
        
        WriteReport();        
    end
    
    FinishLatexDoc();
    dirData = dirResOld;
    recomputeAll = false;      
    
    function InitLatexDoc()
        fileID = fopen(filenameLatexDoc,'w');
    
        fprintf(fileID,'\\documentclass[a4paper,10pt]{article}\n');
        fprintf(fileID,'\\usepackage{listings}\n');
        fprintf(fileID,'\\usepackage[utf8]{inputenc}\n');
        fprintf(fileID,'\\usepackage{graphicx}\n');        
        fprintf(fileID,'\\usepackage[pagebackref,hyperindex=true,colorlinks=true,urlcolor=blue,citecolor=blue,linkcolor=blue]{hyperref}');
        fprintf(fileID,['\\title{ Code for \\lstinline{',dirTest,'} files}\n']);         
        fprintf(fileID,'\\date{\\today}\n');         
        fprintf(fileID,'\\begin{document}\n'); 
        fprintf(fileID,'\\maketitle\n');         
        fprintf(fileID,'\\tableofcontents\n');        
        fclose(fileID);
    end
    
    function AddToLatexDoc(addFile)
         fileID = fopen(filenameLatexDoc,'a');
         fprintf(fileID,['\\section{\\lstinline{',addFile,'}}\n']);      
         
         fprintf(fileID,['{\\scriptsize \\lstinputlisting[language=Matlab]{',replaceBackslash([dirDDFT_2D,'/',addFile]),'}}\n']);
         fclose(fileID);         
    end
    function AddToLatexFig(addFile)
        
         filePath = [dirDDFT_2D_LatexReport filesep addFile];
         fileID = fopen(filenameLatexDoc,'a');                  
         if(exist(filePath,'file'))
             fprintf(fileID,['\\includegraphics[width=8cm]{',replaceBackslash(filePath),'}\n']);
         end
         n = 1;
         while(exist([filePath(1:end-4),num2str(n),'.eps'],'file'))
             fprintf(fileID,['\\includegraphics[width=12cm]{',...
                 replaceBackslash([filePath(1:end-4),num2str(n),'.eps']),'}\n']);
             n = n+1;
         end
         fclose(fileID);         
         
    end
    function FinishLatexDoc()
        fileID = fopen(filenameLatexDoc,'a');
        fprintf(fileID,'\\end{document}\n');
        fclose(fileID); 
        
        %system(['latex ',filenameLatexDoc]);
    end
    function WriteReport()
        tc = 0;
        for j = 1:length(MFiles)
            strfT = MFiles(j).name;
            strfT = strfT(1:end-2);
            if(str2num(Parameters.(strfT)) > 0)
               tc = tc + str2num(Parameters.(strfT));
            end
        end
        
        if(~isfield(Parameters,'DateTime'))
            Parameters.DateTime =  datestr(now);
        end
        Struct2File([dirDDFT_2D,'/TestAll_report.txt'],Parameters,comments);

    end
    function outStr = replaceBackslash(tStr)
        special = '\';
        %tStr = 'Hi, I''m a Big (Not So Big) MATLAB addict; Since my school days!';

        outStr = '';
        for l = tStr
            if (length(find(special == l)) > 0)
                outStr = [outStr, '/'];
            else
                outStr = [outStr, l];
            end
        end
    end                        
end
