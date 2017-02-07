function index = LoadIndexFiles(dirIndexFiles)   

    index         = {};
    TxtFiles      = dir(fullfile(dirIndexFiles,'*.txt'));
    
    [~,dateOrder]  = sort([TxtFiles.datenum],'descend');
    TxtFiles       = TxtFiles(dateOrder);
    
    count = 1;
    
    for i = 1:length(TxtFiles)
        
        st    = File2Struct([dirIndexFiles filesep TxtFiles(i).name]);
        if(isstruct(st))
            index{count} = st;
            count        = count + 1;
        end
    end
    
    
end