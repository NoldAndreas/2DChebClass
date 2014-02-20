function SaveCodeVersion()
    global dirData
    AddPaths();
    time = clock();
	CodeName   = (['2DChebCode_',num2str(time(1)),'_'... % Returns year as character
               num2str(time(2)),'_'...     % Returns month as character
               num2str(time(3)),'_'...     % Returns day as char
               num2str(time(4)),'_'...     % returns hour as char..
               num2str(time(5)),'.zip']);  %returns minute as char                                   
    fullCodePath = [dirData filesep 'CodeVersions' filesep CodeName];
    zip(fullCodePath,pwd);
    disp(['New code version saved in ',fullCodePath])
end