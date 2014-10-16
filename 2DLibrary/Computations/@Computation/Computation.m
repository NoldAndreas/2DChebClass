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
                SaveConfig(this,config,'Configurations');                                                          
                this.optsNum         = config.optsNum;
                this.optsPhys        = config.optsPhys;    
            end
        end
        configName = SaveConfig(this,configuration,configDir)    
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
    end
end