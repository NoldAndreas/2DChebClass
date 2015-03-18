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
        Struct2File(filename,RemovePts(index{ind}));        
        disp(['Pts cleared from found in ',filename,' ...']);
    end
    
    function s = RemovePts(s)
        names = fieldnames(s);
        for k = 1:length(names)            
            if(~isempty(strfind(names{k},'Pts_')))
                s = rmfield(s,[names{k}]);              
            end
        end   
    end
end