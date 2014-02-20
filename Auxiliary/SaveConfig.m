function configName = SaveConfig(configuration,configDir)    
    global dirData

    configDirFull = [dirData filesep configDir];
	if(~exist(configDirFull,'dir'))            
        disp('Folder not found. Creating new path..');            
        mkdir(configDirFull);
    end
   
   configName   = (['Config_',getTimeStr()]);

    
    index         = LoadIndexFiles([dirData filesep configDir]);
    resave        = true;
    config_Flat   = FlattenStructure(configuration,10,'AllStr');
    for i=1:length(index)            
        if(isequalStruct(config_Flat,index{i})) 
            filepath = [configDirFull filesep index{i}.Filename];
            %load(filepath);
            disp(['Using configuration ',filepath,' ...']);
            resave = false;
        end
    end
    if(resave)
        save([configDirFull filesep configName '.mat'],'configuration');
        config_Flat.Filename = [configName,'.mat'];
        
        codename = GetCodeVersionName();
        Struct2File([configDirFull filesep configName '.txt'],config_Flat,...
            ['Configuration ', configName, ' , latest saved Code Version: ',codename]);
        disp(['Configuration saved in ' configDirFull filesep configName '.txt']);
    end

end