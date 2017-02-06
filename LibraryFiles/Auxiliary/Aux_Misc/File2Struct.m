function st = File2Struct(filename)

%    filename = 'Data.txt';

    st     = struct();    
    fileID = fopen(filename);
    
    tline = fgets(fileID);
    line_no = 1;
    while ischar(tline)
                            
        if(tline(1) ~= '#')
            ind        = strfind(tline,':');            
            name       = strtrim(tline(1:ind-1));
            st.(name)  = strtrim(tline(ind+1:end-1));
 %       elseif(line_no == 2)
 %           st.SecondLine = tline(2:end-1);
        end
        
        tline   = fgets(fileID);
        line_no = line_no + 1;
    end
    
    fclose(fileID);

end