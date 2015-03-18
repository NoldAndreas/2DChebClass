function Struct2File(filename,st,comment_str)
    %Saves elements of Structure in File
    %Numbers will be saved up to 14 digits
    
    if((nargin < 3) && isfield(st,'SecondLine'))
        comment_str = st.SecondLine;
        st          = rmfield(st,'SecondLine');
    end

    st = FlattenStructure(st);    
    

    fileID = fopen(filename,'w');
    
    %Print comment:
    fprintf(fileID,'# ************************************************************ \n');
    fprintf(fileID,'# %s \n',comment_str);
    fprintf(fileID,'# ************************************************************ \n');
    
    
    names = fieldnames(st);    
    
    for i = 1:length(names)
        if(isnumeric(st.(names{i})))
            str = num2str(st.(names{i}),12);        
         elseif(iscell(st.(names{i})))
             str = cell2str(st.(names{i}));
        else
            str = st.(names{i});
        end
        fprintf(fileID,'%s : %s\n',names{i},str);
    end
    fclose(fileID);
    
end