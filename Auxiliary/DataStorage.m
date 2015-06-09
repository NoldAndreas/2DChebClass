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
    [h1s,Parameters.Function] = fileparts(func2str(func));
    
	if(isfield(Parameters,'Identifier'))
        identifier = Parameters.Identifier;        
    else
        identifier = [];
    end
    if(isfield(Parameters,'Comments'))
        comments            = Parameters.Comments;            
    else
        comments = '';
    end
	hname               = getFilename(Parameters);        
    filename            = [hname,'.mat'];
    fileParamTxtname    = [hname,'.txt'];            
    %*****************************

    global dirData
    old_dirData = dirData; 
    DataFolder = [dirData filesep nameDir]; 	        
    index      = LoadIndexFiles(DataFolder);   
    
	if(nargin <= 4)
        recompute = [];
    end
    
    if(nargin < 6)
        ignoreList = {};
    end    	
   
    %1st Step: Search for File and load if found
    if((nargin >= 5) && ~(isempty(recompute)) &&  ischar(recompute)) 
        load(recompute);  
        disp(['Stored data from ',recompute,' will be used..']);                        
        recompute = false;
    elseif((nargin >= 5) && ~(isempty(recompute)) && (recompute == 2))        
        %(1)First Option: Select Input File
            [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);            
            load([DataFolder,FileName]);  
            disp(['Stored data from ',[DataFolder,FileName],' will be used..']);                                                
            recompute = false;
    else
        %(2) Second Option: Search file                                
        i          = DataStorage_FindFile(index,Parameters,ignoreList);
        Parameters = index{i};
        
        if(~isempty(i))
            [~,filename,ext]    = fileparts(index{i}.Filename);
            filename            = [filename,ext];
           % filename            = index{i}.Filename;
            Parameters.Filename = [DataFolder filesep filename(1:end-4)];

            fileParamTxtname    = [filename(1:end-4),'.txt'];

            load([DataFolder filesep filename]);
            disp(['Data file found in ',DataFolder filesep filename,' ...']);

            if(recomputeAll)
                recompute = true;
            elseif((isempty(recompute)) && ((~PersonalUserOutput) || (loadAll)))
                recompute = false;
            elseif((isempty(recompute)) && (PersonalUserOutput)) %ask if value of recompute is not given as an input
                no = fprintf('Do you want to load data (press "l"), recompute Data (press any other key), or wait for 2 seconds.');        
                if(getkeywait(2) == 108)
                    [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);            
                    load([DataFolder,FileName]);  
                    disp(['Stored data from ',[DataFolder,FileName],' will be used..']);                        

                    recompute = false;

                elseif(getkeywait(2) == -1)
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
        end     
    
        if(~exist('Data','var'))
            no = fprintf('Do you want to load data (press "l"), or wait for 2 seconds.');        
            if(getkeywait(2) == 108)
                [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);            
                load([DataFolder,FileName]);  
                disp(['Stored data from ',[DataFolder,FileName],' will be used..']);                        
                
                Parameters.Filename =  [DataFolder,FileName];

                recompute = false;
            else
                for ih = 1:no
                    fprintf('\b');
                end  
                recompute = true;       
            end
        end         
        
        %2b) If File not found: Compute Data and Save Data
        %if(~exist('Data','var') || ((nargin == 5) && recompute == true))
        if(recompute)
           startComp = tic;
           fprintf(['Computing data, starting at ', datestr(now), '...\n']);
           if(nargout(func)>=2)
                [Data,Parameters.Results] = func(OriginalParameters,OtherInputs);
           else
                Data     = func(OriginalParameters,OtherInputs);
           end
           t                               = toc(startComp);                
           
           if(~exist(DataFolder,'dir'))            
                disp('Folder not found. Creating new path..');            
                mkdir(DataFolder);
           end
           save([DataFolder filesep filename],'Data');
           disp(['Data saved in ',[DataFolder filesep filename],'.']);            
            
           Parameters.Results.comptime_sec = t;
           Parameters.Results.comp_at      = datestr(now);           
           if(isnumeric(Data) && isscalar(Data) && (Data == 0))
                disp(['Error after ',sec2hms(t), ' (hrs:min:sec) ']);
                comments = ['Computation failed! ',comments];
           else
                disp(['Data recomputed: ',sec2hms(t),' (hrs:min:sec)']);
           end
           Parameters.Results.comments  = comments; 
           Parameters.Filename          = filename;          
            
           Struct2File([DataFolder filesep fileParamTxtname],Parameters,...
                ['Computed at: ',datestr(now),...
                ' Computation time: ',sec2hms(t), ' (hrs:min:sec) ',comments]);
            
            %disp(['Index file saved in ',[DataFolder filesep fileParamTxtname],'..']);            
            Parameters.Filename = [DataFolder filesep Parameters.Filename(1:end-4)];
         end       
    end
    ChangeDirData(old_dirData);

    
    function name = getFilename(params)        
        if(isempty(identifier))
            name   = '';
        elseif(isfloat(params.(identifier)))
            name = [identifier,'_',num2str(params.(identifier)),'_'];
        else
            name = [identifier,'_',params.(identifier),'_'];
        end
        name   = ([name,getTimeStr()]);        
    end

end