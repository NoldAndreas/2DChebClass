function [Data,recompute,Parameters] = DataStorage(nameDir,func,Parameters,OtherInputs,recompute,ignoreList)
%Input: 
%Name of Folder,List of Parameters
%if file is to be written: Data
%
%!!! Attention: As Parameters, numbers will be trimmed to 12 digits !!!

    %%*****************************
    %Check if there is a personal user!
    global PersonalUserOutput recomputeAll loadAll
    %*****************************
    %Trim numbers of parameters
    OriginalParameters = Parameters;
    Parameters         = FlattenStructure(Parameters,10,'AllStr');    
    
	if(isfield(Parameters,'Identifier'))
        identifier = Parameters.Identifier;
        Parameters = rmfield(Parameters,'Identifier');
    else
        identifier = [];
    end
    %*****************************

    global dirData
    DataFolder = [dirData filesep nameDir]; 	        
    index      = LoadIndexFiles(DataFolder);   
    
	if(nargin == 4)
        recompute = [];
    end
   
   
    %1st Step: Search for File and load if found
    if((nargin >= 5) && ~(isempty(recompute)) && (recompute == 2))        
        %(1)First Option: Select Input File
            [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);            
            load([DataFolder,FileName]);  
            disp(['Stored data from ',[DataFolder,FileName],' will be used..']);                        
                        
%            for i=1:length(index)
%                [~,name,ext] = fileparts(index(i).Filename);
%               if(strcmp(FileName,[name,ext]))
%                    Parameters = index(i).Parameters;                    
%                    break;
%                end
%            end            
            recompute = false;
    else
        %(2) Second Option: Search file
        if(isfield(Parameters,'Comments'))
            comments            = Parameters.Comments;
            Parameters          = rmfield(Parameters,'Comments');
        else
            comments = '';
        end
                        
        Parameters.Function = func2str(func);
        hname               = getFilename(Parameters);        
        filename            = [hname,'.mat'];
        fileParamTxtname    = [hname,'.txt'];
            
        [h1s,Parameters.Function] = fileparts(Parameters.Function);
        if(nargin >= 6)
            Parameters_comp = RemoveIgnoreFromStruct(Parameters,ignoreList);
        else
            Parameters_comp = Parameters;
        end
        for i=1:length(index)
            if(isfield(index{i},'Filename'))
                par_i = rmfield(index{i},'Filename');
            else
                par_i = index{i};
            end
                        
            [h1s,par_i.Function]      = fileparts(par_i.Function);
              
%             comp_struct(Parameters,par_i)
%             pause
            if(nargin >= 6)
                par_i_comp = RemoveIgnoreFromStruct(par_i,ignoreList);
            else
                par_i_comp = par_i;
            end

            if(isequaln(Parameters_comp,par_i_comp))%index{i}.Parameters))
                
                filename            = index{i}.Filename;
                Parameters.Filename = filename(1:end-4);
                fileParamTxtname    = [filename(1:end-4),'.txt'];
                
                load([DataFolder filesep filename]);
                disp(['Data file found in ',DataFolder filesep filename,' ...']);

                if(recomputeAll)
                    recompute = true;
                elseif((isempty(recompute)) && ((~PersonalUserOutput) || (loadAll)))
                    recompute = false;
                elseif((isempty(recompute)) && (PersonalUserOutput)) %ask if value of recompute is not given as an input
                    no = fprintf('Do you want to recompute Data (press any key), or wait for 2 seconds.');        
                    if(getkeywait(2) == -1)
                        for ih = 1:no
                            fprintf('\b');
                        end                        
                        recompute = false;
                    else
                        for ih = 1:no
                            fprintf('\b');
                        end
                        no = fprintf('Thanks. Data will be recomputed.\n');
                        recompute = true;
                    end                
                end
                if(~recompute)
                    %append configuration to index file
                     fileID = fopen([DataFolder filesep fileParamTxtname],'a');
                     fprintf(fileID,'# %s \n',['loaded on ' datestr(clock) ' : ' comments]);
                     fclose(fileID);
                end
                break;
            end
        end     
    
        if(~exist('Data','var'))                           
            recompute = true;       
        end         
        
        %2b) If File not found: Compute Data and Save Data
        %if(~exist('Data','var') || ((nargin == 5) && recompute == true))
        if(recompute)
           startComp = tic;
           fprintf(['Computing data, starting at ', datestr(now), '...\n']);
           Data     = func(OriginalParameters,OtherInputs);
           t        = toc(startComp);     
           if(isnumeric(Data) && isscalar(Data) && (Data == 0))
                disp(['Error after ',sec2hms(t), ' (hrs:min:sec) ']);
                comments = ['Computation failed! ',comments];
           else
                disp(['Data recomputed: ',sec2hms(t),' (hrs:min:sec)']);
           end

           Parameters.Filename = filename;
           
           if(~exist(DataFolder,'dir'))            
                disp('Folder not found. Creating new path..');            
                mkdir(DataFolder);
           end
            save([DataFolder filesep filename],'Data');
            
            Struct2File([DataFolder filesep fileParamTxtname],Parameters,...
                ['Computed at: ',datestr(now),...
                ' Computation time: ',sec2hms(t), ' (hrs:min:sec) ',comments]);

            disp(['Data saved in ',[DataFolder filesep filename],'.']);            
            %disp(['Index file saved in ',[DataFolder filesep fileParamTxtname],'..']);            
            Parameters.Filename = Parameters.Filename(1:end-4);
         end       
    end

    
    function name = getFilename(params)
        time   = clock;
%        fnames = fieldnames(params);       
%        name   = '';%[func2str(func),'_'];
%         if(isempty(identifier))
%             k = 0; kMax = 3;
%             for j = 1:length(fnames)            
%                 if(isfloat(params.(fnames{j})))
%                     s    = num2str(params.(fnames{j}));
%                     if(length(s) < 6)
%                         name = [name,fnames{j},'_',s(:)','_'];
%                         k = k +1;
%                     end
%                 end
%                 if((k > kMax) || length(name) > 35)
%                     break;
%                 end
%             end
%         else
        if(isempty(identifier))
            name   = '';
        elseif(isfloat(params.(identifier)))
            name = [identifier,'_',num2str(params.(identifier)),'_'];
        else
            name = [identifier,'_',params.(identifier),'_'];
        end
%       end
        name   = ([name,getTimeStr()]);        
    end

%     function eq = isequalStruct(s1,s2)
%         
%         eq       = true;        
%         f_names1 = fieldnames(s1);
%         tol      = 1e-12;
%         
%         for j = 1:length(f_names1)
%             fld = f_names1{j};
%             if(~isfield(s2,fld))
%               eq = false; return;
%             else
%                 if(isnumeric(s1.(fld)))
%                     if(abs(s1.(fld)- s2.(fld)) > tol)
%                         eq = false; return;
%                     end
%                 else
%                     if(~strcmp(s1.(fld),s2.(fld)))
%                         eq = false; return;
%                     end
%                 end
%             end            
%         end
%        
%     end

    function s = RemoveIgnoreFromStruct(s,IgnoreList)
        for j = 1:length(IgnoreList)
            for k = 1:length(s)
                
            end
            
            if(isfield(s,IgnoreList{j}))
                s = rmfield(s,IgnoreList{j});
            
            end
        end
    end

end


%             %Save text file with Parameters
%             [names,values] = fn_structdispList(newEntry); 
%             
%             fileID = fopen(fileParamTxtname,'w');       
%             fprintf(fileID,'****** Data File *******\n');
%             fprintf(fileID,'Filename: %s\n',filename);
%             fprintf(fileID,'************************\n');
%             for i = 1:length(names)      
%                 if(ischar(values{i}))
%                     fprintf(fileID,'%s : %s\n',names{i},values{i});
%                 else
%                     if(length(values{i}) == 1)
%                         fprintf(fileID,'%s : %f\n',[names{i},values{i}]);                
%                     else
%                         fprintf(fileID,'%s : ',[names{i}]);                
%                         for l = 1:length(values{i})
%                             fprintf(fileID,'%f , ',values{i}(l));
%                         end
%                         fprintf(fileID,'\n');  
%                     end
%                 end
%             end
%             fclose(fileID);