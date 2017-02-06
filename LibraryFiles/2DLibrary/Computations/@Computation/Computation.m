classdef Computation < handle
    properties (Access = public)            
        optsNum,optsPhys
        configName
        IDC
        
        %subArea-parameters and function
        subArea,doSubArea,Int_of_path,IP
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
        
        function PreprocessIDC(this)
            optsNum  = this.optsNum;
            shape    = optsNum.PhysArea;            
            %*****************************************************
            %Special Cases - correct configuration, 
            %   account for outdated configuration files
            if(isfield(optsNum,'V2Num'))
                shape.Conv = optsNum.V2Num;
            else                
                shape.Conv = [];
            end
            
            % Special Case: HalfSpace_FMT
            if(strcmp(optsNum.PhysArea.shape,'HalfSpace_FMT') || ...
               strcmp(optsNum.PhysArea.shape,'InfCapillary_FMT'))
                shape.R = this.optsPhys.sigmaS/2;
            end
            %*****************************************************
            
            % Construct main object for geometry IDC
            shapeClass = str2func(optsNum.PhysArea.shape);
            this.IDC   = shapeClass(shape);
            if(isfield(optsNum,'PlotArea'))
                this.IDC.ComputeAll(optsNum.PlotArea);             
            elseif(isfield(optsNum,'PlotAreaCart'))
                this.IDC.ComputeAll();
                this.IDC.InterpolationPlotCart(optsNum.PlotAreaCart,true);
                this.IDC.InterpolationPlotFlux(optsNum.PlotAreaCart);
            else
                this.IDC.ComputeAll();
            end
        end
        function Preprocess(this)                                                 

            if(~isfield(this.optsPhys,'sigmaS') && isfield(this.optsPhys,'V2') && isfield(this.optsPhys.V2,'sigmaS'))
                this.optsPhys.sigmaS = this.optsPhys.V2.sigmaS;
            end                            
            
            PreprocessIDC(this);
            Preprocess_SubArea(this);
        end        
        function Preprocess_SubArea(this)
            
            optsNum  = this.optsNum;
            optsPhys = this.optsPhys;

            this.doSubArea = isfield(optsNum,'SubArea');

            if(this.doSubArea)    
                subshapeClass  = str2func(optsNum.SubArea.shape);
                this.subArea   = subshapeClass(optsNum.SubArea);
                
                plotSubShape   = optsNum.SubArea;
                plotSubShape.N = [80,80];
                this.subArea.ComputeAll(plotSubShape); 
                
                this.IP            = this.IDC.SubShapePtsCart(this.subArea.GetCartPts());
                this.Int_of_path   = this.subArea.IntFluxThroughDomain(100)*blkdiag(this.IP,this.IP);
            else
                this.Int_of_path   =  zeros(1,2*this.IDC.M);
            end
        end        

        function conf = GetConfig(this)
            conf = struct('optsNum',this.optsNum,'optsPhys',this.optsPhys);
        end
        function configName = SaveConfig(this,configuration)    
        
            global dirDataOrg
            if(nargin < 2)
                optsNum  = this.optsNum;
                optsPhys = this.optsPhys;
                configuration = v2struct(optsNum,optsPhys);
            end
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
        function fullName = SaveCurrentFigure(this,filename)                                    
            s.optsNum     = this.optsNum;
            s.optsPhys    = this.optsPhys;
            s.configName  = this.configName;
            fullName      =  SaveFigure(filename,s);
        end
%         
%         function LoadSnapshot(this)
%             global dirDataOrg
%             
%             datafileName  = [this.configName '_Data'];
%             configDirFull = [dirDataOrg filesep 'Configurations' filesep  datafileName '.mat'];            
%             if(exist(configDirFull,'file'))
%                 load(configDirFull);            
%                 disp(['Data file ',configDirFull,' loaded.']);
%             end
%             
%         end
%         function SaveSnapshot(this)
%             global dirDataOrg
%             
%             datafileName  = [this.configName '_Data'];
%             configDirFull = [dirDataOrg filesep 'Configurations' filesep  datafileName '.mat'];            
%             
%             save(configDirFull,'this');
%             disp(['Data saved in ' configDirFull]);
%         end
        
    end
end