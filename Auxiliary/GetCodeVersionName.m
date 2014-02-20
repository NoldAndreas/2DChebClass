function nameCurrentCodeVersion = GetCodeVersionName()
    global CodeVersionsDir
    ZipFiles      = dir(fullfile([CodeVersionsDir],'*.zip'));
    
    dates = zeros(length(ZipFiles),1);
    for i = 1:length(ZipFiles)
        dates(i) = ZipFiles(i).datenum;
    end
    
    [h1,i] = max(dates);
    nameCurrentCodeVersion = ZipFiles(i).name;
end