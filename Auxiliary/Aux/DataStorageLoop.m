function [Data,recompute,Parameters] = DataStorageLoop(nameDir,func,Parameters,OtherInputs,recompute,noLoops)
%Input: 
%Name of Folder,List of Parameters
%if file is to be written: Data
%
%!!! Attention: As Parameters, numbers will be trimmed to 12 digits !!!

    %*****************************
    %Trim numbers of parameters
    OriginalParameters = Parameters;
    Parameters         = FlattenStructure(Parameters,10,'AllStr');    
%    Parameters         = FlattenStructure(Parameters,10);    
    
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
        %set default value
        recompute = true;
    end
    
    %1st Step: Search for File and load if found
    if((nargin == 5) && (recompute == 2))        
        %(1)First Option: Select Input File
            [FileName,DataFolder] = uigetfile('*.mat',['Select Data File for ',func2str(func)]);            
            load([DataFolder,FileName]);  
            disp(['Stored data from ',[DataFolder,FileName],' will be used..']);                        
                        
            for i=1:length(index)
                [~,name,ext] = fileparts(index(i).Filename);
                if(strcmp(FileName,[name,ext]))
                    Parameters = index(i).Parameters;                    
                    break;
                end
            end            
            recompute = false;
    else
        %(2) Second Option: Search file
        Parameters.Function = func2str(func);
        hname               = getFilename(Parameters);        
        filename            = [hname,'.mat'];
        fileParamTxtname    = [hname,'.txt'];
                
        for i=1:length(index)
            par_i = rmfield(index{i},'Filename');
            if(isequalStruct(Parameters,par_i))%index{i}.Parameters))
                
                filename         = index{i}.Filename;
                fileParamTxtname = [filename(1:end-4),'.txt'];
                
                load([DataFolder filesep filename]);
                disp(['Loading file ',DataFolder filesep filename,' ...']);
                
                if((nargin==4) || isempty(recompute)) %ask if value of recompute is not given as an input
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
        
        %2b) If File not found: Compute Data and Save Data
        %if(~exist('Data','var') || ((nargin == 5) && recompute == true))
        if(recompute)

           if(~exist(DataFolder,'dir'))            
                disp('Folder not found. Creating new path..');            
                mkdir(DataFolder);
           end
           Parameters.Filename = filename;                                 
 
           tic
           disp('Computing data ...')                      
           h = waitbar(0,'Computing Iterations...');
           t = 0;
           for iloop = 1:noLoops
                OriginalParameters.i = iloop;
                [Data(iloop,:),stop] = func(OriginalParameters,OtherInputs);
                OtherInputs = Data(iloop,:)';
                
                save([DataFolder filesep filename],'Data');
                %disp(['Data file saved in ',[DataFolder filesep filename],'..']);                            
                %if(iloop == 1)
                 t  = t + toc;
                 tic
                 Struct2File([DataFolder filesep fileParamTxtname],Parameters,['Computed at: ',datestr(now),' Computation time for ',num2str(iloop),' steps: ',num2str(t),' sec']);                    
                 %disp(['Index file saved in ',[DataFolder filesep fileParamTxtname],'..']);            
                %end
                waitbar(iloop/noLoops,h);
                if(stop)
                    break;
                end
           end
           t = toc;     
           close(h);       
           disp(['Complete Data recomputed: ',num2str(t),' sec']);
           
%            disp(['Data file saved in ',[DataFolder filesep filename],'..']);                                    
%            disp(['Index file saved in ',[DataFolder filesep fileParamTxtname],'..']);            

        end       
    end

    
    function name = getFilename(params)
        time   = clock;
        
        if(isempty(identifier))
            name   = '';
        elseif(isfloat(params.(identifier)))
            name = [identifier,'_',num2str(params.(identifier)),'_'];
        else
            name = [identifier,'_',params.(identifier),'_'];
        end
        
        name   = ([name,getTimeStr()]);        
    end

    function eq = isequalStruct(s1,s2)
        
        eq       = true;        
        f_names1 = fieldnames(s1);
        tol      = 1e-12;
        
        for j = 1:length(f_names1)
            fld = f_names1{j};
            if(~isfield(s2,fld))
              eq = false; return;
            else
                if(isnumeric(s1.(fld)))
                    if(abs(s1.(fld)- s2.(fld)) > tol)
                        eq = false; return;
                    end
                else
                    if(~strcmp(s1.(fld),s2.(fld)))
                        eq = false; return;
                    end
                end
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