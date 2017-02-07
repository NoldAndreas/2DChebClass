function ind = DataStorage_FindFile(index,Parameters,ignoreList)
         
    ignoreList      = [ignoreList,{'Results','Comments'}];                               
    ignoreListSingle = {'configName','Filename','Identifier','SecondLine'};
                               
                               
    Parameters_comp = RemoveIgnoreFromStruct(Parameters,ignoreList);
    
    ind = [];
    for i=1:length(index)        
        par_i = index{i};        

        if(~isfield(par_i,'Function'))
            continue;
        end
        [h1s,par_i.Function]      = fileparts(par_i.Function);

        par_i_comp = RemoveIgnoreFromStruct(par_i,ignoreList);
        
        if(isequaln(Parameters_comp,par_i_comp))%index{i}.Parameters))
            ind = i;
            break;
        end
    end
    
    function s = RemoveIgnoreFromStruct(s,IgnoreList)
        
        %1st step: remove the fields of s that appear in IgnoreList
        jj = isfield(s,ignoreListSingle);
        s  = rmfield(s,ignoreListSingle(jj));
        
        jj = isfield(s,IgnoreList);
        s  = rmfield(s,IgnoreList(jj));
        IgnoreList = IgnoreList(~jj);
        
        names = fieldnames(s);
        for k = 1:length(names)
            for j = 1:length(IgnoreList)            
               if(~isempty(strfind(names{k},[IgnoreList{j},'_'])))
                   %s.(names{k}) = [];
                   s = rmfield(s,[names{k}]);              
               end
            end           
        end      
        
    end
end