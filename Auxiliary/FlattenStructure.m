function C = FlattenStructure(A,decimalCut,AllStr)
% 
%(1) Flattens a nested, recursive Structure
%(2) Transforms all vector entries to scalar ones
%(3) Cuts down numbers to decumalCut decimal places
%
    if(nargin == 2)
        AllStr = '';
    end
        
    C = struct();
    names = fieldnames(A);

    for i=1:length(names)  
        if(~isstruct(A.(names{i})))
            C.(names{i}) = A.(names{i});
        else
           Ctemp = FlattenStructure(A.(names{i}));
           C     = concatenateStructuresChangeName(C,Ctemp,names{i});% [C,Ctemp];       
        end

    end
    
    %(2) Transform vectors and cells to multiple entries
    f_names     = fieldnames(C);
    for i=1:length(f_names)
        val  = C.(f_names{i});        
        if(isnumeric(val))
            [s1,s2] = size(val);
            if((s1>1) || (s2 > 1))
                for j1 = 1:s1
                    for j2 = 1:s2                        
                        if((s1 > 1) && (s2 > 1))
                            newname = [f_names{i},'_',num2str(j1),'_',num2str(j2)];
                        elseif((s1 == 1) || (s2 == 1))
                            newname = [f_names{i},'_',num2str(max(j1,j2))];                           
                        end                        
                        C.(newname) = val(j1,j2);
                    end
                end
                C = rmfield(C,f_names{i});
            end
        elseif(iscell(val))
            [s1,s2] = size(val);
            if(s1==1 && s2==1)
                C.(f_names{i}) = C.(f_names{i}){1};
            else
                for j1 = 1:s1
                    for j2 = 1:s2                        
                        if((s1 > 1) && (s2 > 1))
                            newname = [f_names{i},'_',num2str(j1),'_',num2str(j2)];
                        elseif((s1 == 1) || (s2 == 1))
                            newname = [f_names{i},'_',num2str(max(j1,j2))];                           
                        end                        
                        C.(newname) = val{j1,j2};
                    end
                end
                C = rmfield(C,f_names{i});
            end
        end
    end
     
    %(3) Also cut down numbers to decimalCut decimal places    
    if(nargin > 1)        
        f_names     = fieldnames(C);
        for i=1:length(f_names)
            % deal with logical values
            if(islogical(C.(f_names{i})))
                C.(f_names{i}) = double(C.(f_names{i}));
            end
            if(isempty(C.(f_names{i})))
                C.(f_names{i}) = '[]';
            elseif(isnumeric(C.(f_names{i})))
                C.(f_names{i}) = roundsd(C.(f_names{i}),decimalCut);
            end
        end
    end
    
    if((nargin > 2) && strcmp(AllStr,'AllStr'))
        for i=1:length(f_names)
            if(isnumeric(C.(f_names{i})))
                C.(f_names{i}) = num2str(C.(f_names{i}));
            end
        end
    end
   
    function snew = concatenateStructures(s1,s2)
        snew = struct();
        
        names1 = fieldnames(s1);
        for i1 = 1:length(names1)  
            snew.(names1{i1}) = s1.(names1{i1});
        end
        
        names2 = fieldnames(s2);
        for i2 = 1:length(names2)  
            snew.(names2{i2}) = s2.(names2{i2});
        end        
    
    end

    function snew = concatenateStructuresChangeName(s1,s2,str2)
        snew = struct();
        
        names1 = fieldnames(s1);
        for i1 = 1:length(names1)  
            snew.(names1{i1}) = s1.(names1{i1});
        end
        
        names2 = fieldnames(s2);
        for i2 = 1:length(names2)  
            snew.([str2,'_',names2{i2}]) = s2.(names2{i2});
        end        
    
    end

end
% 
% function s = FlattenStructure(sIn)
% %Flattens a recursive structure
%     
%     names = fieldnames(sIn);
%     for i = 1:length(names)
%                 
%         if(isstruct(sIn.(names{i})))
%             disp(sIn.(names{i}));
%         end
%         
%     end        
%     
%     function smod = ModifyStructureNames(sIn)
%         
%     end
% 
% end