function ClearPts(nameDir,func,Parameters,ignoreList)  
    if(~isfield(Parameters,'Pts'))
        disp('Pts cannot be cleared as not part of input structure');
        return;
    end
    
    if(~isfield(Parameters,'physArea'))
        disp('Pts cannot be replaced as physArea not part of input structure');
        return;
    end               
    
    global dirData
    DataFolder = [dirData filesep nameDir];
    index      = LoadIndexFiles(DataFolder);
    
    [h1s,Parameters.Function] = fileparts(func2str(func));
    Parameters         = FlattenStructure(Parameters,10,'AllStr');    
    ind                = DataStorage_FindFile(index,Parameters,ignoreList);
    
    if(~isempty(ind))
        filename         = [DataFolder filesep index{ind}.Filename(1:end-4),'.txt'];
        s                = RemovePts(index{ind});
        s                = AddPhysArea(s,Parameters);
        Struct2File(filename,s);        
        disp(['Pts cleared from ',filename,' ...']);
    end
    
    
    function sOut = RemovePts(s)
        names = fieldnames(s);
        for k = 1:length(names)            
            if(isempty(strfind(names{k},'Pts_')))
                sOut.(names{k}) = s.(names{k});                
            end
        end   
    end        

    function s = AddPhysArea(s,parameters)
                
        names = fieldnames(s);
        for k = 1:length(names)            
            if(~isempty(strfind(names{k},'physArea')))
                return;
            end
        end   
        
        %if no physArea found: add
        names = fieldnames(parameters);
        for k = 1:length(names)            
            if(~isempty(strfind(names{k},'physArea')))
                s.(names{k}) = parameters.(names{k});
            end
        end   
    end
        
end