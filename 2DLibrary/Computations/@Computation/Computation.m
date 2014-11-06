classdef Computation < handle
    properties (Access = public)            
        optsNum,optsPhys
        configName
    end
    methods (Access = public)   
        function this = Computation(config)            
            if(isempty(config))
                LoadConfig(this);
            else
                SaveConfig(this,config);                                                          
                this.optsNum         = config.optsNum;
                this.optsPhys        = config.optsPhys;    
            end
        end
        function configName = SaveConfig(this,configuration)    
        
            global dirDataOrg
            configDir = 'Configurations';

            configDirFull = [dirDataOrg filesep configDir];
            if(~exist(configDirFull,'dir'))            
                disp('Folder not found. Creating new path..');            
                mkdir(configDirFull);
            end
                        
            index         = LoadIndexFiles([dirDataOrg filesep configDir]);
            resave        = true;
            config_Flat   = FlattenStructure(configuration,10,'AllStr');
            for i=1:length(index)            
                if(isequalStruct(config_Flat,index{i})) 
                    filepath = [configDirFull filesep index{i}.Filename];
                    %load(filepath);
                    disp(['Using configuration ',filepath,' ...']);
                    this.configName = index{i}.Filename(1:end-4);
                    resave = false;
                    LoadSnapshot(this);
                end
            end
            if(resave)
                configName   = (['Config_',getTimeStr()]);
                save([configDirFull filesep configName '.mat'],'configuration');
                config_Flat.Filename = [configName,'.mat'];

                Struct2File([configDirFull filesep configName '.txt'],config_Flat,...
                    ['Configuration ', configName]);
                disp(['Configuration saved in ' configDirFull filesep configName '.txt']);
                this.configName = configName;
            end
        end
        function LoadConfig(this)
            global dirData

            [configIn,DataFolder] = uigetfile([dirData filesep 'Configurations' filesep '*.mat'],['Select Config File']);
            load([DataFolder,configIn]);
            disp(['Loading configuration ',[DataFolder,configIn],' ...']);
            this.configName = configIn(1:end-4); 
            this.optsNum    = configuration.optsNum;
            this.optsPhys   = configuration.optsPhys;
        end
        function SaveCurrentFigure(this,filename)
            global dirData
            print2eps([dirData filesep filename],gcf);
            saveas(gcf,[dirData filesep filename '.fig']);            
            disp(['Figures saved in ',dirData filesep filename '.fig/eps']);
        end
        
        function LoadSnapshot(this)
            global dirDataOrg
            
            datafileName  = [this.configName '_Data'];
            configDirFull = [dirDataOrg filesep 'Configurations' filesep  datafileName '.mat'];            
            if(exist(configDirFull,'file'))
                load(configDirFull);            
                disp(['Data file ',configDirFull,' loaded.']);
            end
            
        end
        function SaveSnapshot(this)
            global dirDataOrg
            
            datafileName  = [this.configName '_Data'];
            configDirFull = [dirDataOrg filesep 'Configurations' filesep  datafileName '.mat'];            
            
            save(configDirFull,'this');
            disp(['Data saved in ' configDirFull]);
        end
        
    end
end