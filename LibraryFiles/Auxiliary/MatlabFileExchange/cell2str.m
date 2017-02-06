function string = cell2str(cellstr)

    string = '';
    for i = 1:length(cellstr)        
        if iscell(cellstr{i})
            str = cell2str(cellstr{i});
        elseif ischar(cellstr{i})
            str = cellstr{i};
        end    
        string = [string str '_'];
    end
    
    string = string(1:end-1);
end