function st = File2Struct(filename)

%    filename = 'Data.txt';

    st     = struct();    
    fileID = fopen(filename);
    
    tline = fgets(fileID);
    while ischar(tline)
                            
        if(tline(1) ~= '#')
            [name,remain]  = strtok(tline);
            [h1,remain]     = strtok(remain);
            data           = strtok(remain);

%            if(isempty(str2num(data)))
            st.(name) = data;            
%            else
%                st.(name) = str2num(data);            
%            end
        end
        
        tline = fgets(fileID);
    end
    
    fclose(fileID);

end