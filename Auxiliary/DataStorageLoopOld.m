function [Data,recompute] = DataStorageLoop(nameDir,func,Parameters,OtherInputs,recompute,noLoops)
%Input: 
%Name of Folder,List of Parameters
%if file is to be written: Data

    %*****************************
    %Trim numbers of parameters
    OriginalParameters = Parameters;
    Parameters         = FlattenStructure(Parameters,10);    
    
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
        
    Parameters.Function = func2str(func);   
	hname               = getFilename(Parameters);        
    filename            = [hname,'.mat'];
    fileParamTxtname    = [hname,'.txt'];
    
    %1st Step: Search for File and load if found
    if((nargin == 5) && (recompute == 2))
            [FileName,PathName] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);
            load([PathName,FileName]);            
            disp(['Stored data from ',[PathName,FileName],' will be used..']);
            recompute = false;
    else
           for i=1:length(index)
            par_i = rmfield(index{i},'Filename');
            if(isequalStruct(Parameters,par_i))%index{i}.Parameters))
                
                filename         = index{i}.Filename;
                fileParamTxtname = [filename(1:end-4),'.txt'];
                
                load([DataFolder filesep filename]);
                disp(['Loading file ',DataFolder filesep filename,' ...']);
                
                if(nargin==4) %ask if value of recompute is not given as an input
                    disp('Do you want to recompute the Data (press any key), otherwise wait for 2 seconds.');        
                    if(getkeywait(2) == -1)
                        disp('Stored data will be used');
                        recompute = false;
                    else
                        disp('Thanks. Data will be recomputed.');
                        recompute = true;
                    end                
                end
                break;
            end
        end     
    
        if(~exist('Data','var'))               
            disp('Recomputing data...');
            recompute = true;       
        end 
    end
    
    %2b) If File not found: Compute Data and Save Data
    %if(~exist('Data','var') || ((nargin == 5) && recompute == true))
    if(recompute)
        
       newEntry = struct('Parameters',Parameters,...
                         'Filename',filename);
       index    = [index;newEntry];   
       
       fpathstr = fileparts(filename);
       if(~exist(fpathstr,'dir'))            
            disp('Folder not found. Creating new path..');            
            mkdir(fpathstr);
       end       
        
       tic            
       h = waitbar(0,'Computing Iterations...');
       for iloop = 1:noLoops
            Parameters.i = iloop;
            [Data(iloop,:),stop] = func(OriginalParameters,OtherInputs);
            save(filename,'Data');
            %disp('Data file saved..');
            if(iloop == 1)
                save(IndexFile,'index'); 
                disp('Index file saved..');
            end
            waitbar(iloop/noLoops,h);  
            if(stop)
                break;
            end
       end
       t        = toc;     
       close(h);       
       disp(['Complete Data recomputed: ',num2str(t),' sec']);                       

    end       

    
    function name = getFilename(params)
        time   = clock;
        fnames = fieldnames(params);
        
        name   = [func2str(func),'_'];
        k = 0; kMax = 3;
        for j = 1:length(fnames)            
            if(isfloat(params.(fnames{j})))
                s    = num2str(params.(fnames{j}));
                if(length(s) < 6)
                    name = [name,fnames{j},'_',s(:)','_'];
                    k = k +1;
                end
            end
            if((k > kMax) || length(name) > 35)
                break;
            end
        end
        
        name   = ([name,'_',num2str(time(1)),'_'... % Returns year as character
                    num2str(time(2)),'_'... % Returns month as character
                    num2str(time(3)),'_'... % Returns day as char
                    num2str(time(4)),'_'... % returns hour as char..
                    num2str(time(5)),... %returns minute as char                    
                    '.mat']);        
    end
    

end