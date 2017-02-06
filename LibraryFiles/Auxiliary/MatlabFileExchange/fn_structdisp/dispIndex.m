function dispIndex(index)
    for j = 1:length(index)
        [names(:,j),values(:,j)] = fn_structdispList(index(j).Parameters);        
    end
    
    headline = {'Parameter'};
    for j = 1:length(index)
        [~, fname, ~] = fileparts(index(j).Filename);
        headline{j+1} = fname;
    end
    
    disp(headline);
    disp([headline;names(:,1),values]);
    
    
end