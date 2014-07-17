classdef DDFT_2D < handle
    
    properties (Access = public)            
        optsNum,optsPhys        
        
        IDC,subArea
        IntMatrV2  % integration matrix for mean field two-particle interactions
        IntMatrHI  % integration matrices for hydrodynamic interactions
        IntMatrFex % integration matrix for FMT (hard-sphere) interactions
        Int_of_path
        IP
        
        Vext,Vext_grad,VAdd
        
        x_eq,mu
        dynamicsResult
        
        doHI,doSubArea,do2Phase
    end
    
    methods (Access = public)          
        function this = DDFT_2D(configuration)
            
            this.optsNum         = configuration.optsNum;
            this.optsPhys        = configuration.optsPhys;    
            
            if(~isfield(this.optsPhys,'D0') && isfield(this.optsPhys,'mS') && isfield(this.optsPhys,'gammaS'))
                this.optsPhys.D0 = 1./(this.optsPhys.gammaS.*this.optsPhys.mS);
            elseif(~isfield(this.optsPhys,'D0'))
                this.optsPhys.D0 = 1;
            end                
            
        end        
        function Preprocess(this)
            
            optsNum  = this.optsNum;
            optsPhys = this.optsPhys;
            
            shape    = optsNum.PhysArea;
            if(isfield(optsNum,'V2Num'))
                shape.Conv = optsNum.V2Num;
            else                
                shape.Conv = [];
            end
            
            % Special Case: HalfSpace_FMT
            if(strcmp(optsNum.PhysArea.shape,'HalfSpace_FMT'))
                shape.R = this.optsPhys.sigmaS/2;
            end
            
            % Construct main object for geometry IDC
            shapeClass = str2func(optsNum.PhysArea.shape);
            this.IDC   = shapeClass(shape);            
            this.IDC.ComputeAll(optsNum.PlotArea);             
            
            % Determine hard-sphere contribution to bulk free energy
            if(~isfield(this.optsPhys,'HSBulk'))
                if(isfield(this.optsPhys,'FexNum'))
                    this.optsPhys.HSBulk = ['FexBulk_',this.optsPhys.FexNum.Fex];
                else
                    this.optsPhys.HSBulk = 'Fex_ZeroMap';
                end
            end
            
            % Determine saturation point if a 2-Phase system is given
            if(isfield(this.optsPhys,'V2') && ~strcmp(this.optsPhys.HSBulk,'Fex_ZeroMap'))
                [this.optsPhys.rhoGas_sat,...
                     this.optsPhys.rhoLiq_sat,...
                    this.optsPhys.mu_sat,p] = BulkSatValues(this.optsPhys);
                GetCriticalPoint(this.optsPhys);
                % BulkPhaseDiagram(this.optsPhys);
                this.do2Phase = true;
            end

            % Determine number of species
            if(isfield(optsPhys,'nParticlesS'))
                this.optsPhys.nSpecies = length(optsPhys.nParticlesS); 
            elseif(isfield(optsPhys,'Dmu'))
                this.optsPhys.nSpecies = length(optsPhys.Dmu);
            elseif(isfield(optsPhys,'mu'))
                this.optsPhys.nSpecies = length(optsPhys.mu);
            end
                                                             
            Preprocess_HardSphereContribution(this);
            Preprocess_MeanfieldContribution(this);
            Preprocess_HIContribution(this);           
            Preprocess_ExternalPotential(this);                        
            Preprocess_SubArea(this);
            
        end        
        
        function y0 = getInitialGuess(this,rho_ig)
            
          optsNum  = this.optsNum;
          optsPhys = this.optsPhys;          
          
          if(nargin >= 2)
                if(isscalar(rho_ig))
                    rho_ig = rho_ig*ones(this.IDC.M,1);
                end
                y0 = this.optsPhys.kBT*log(rho_ig) + this.Vext;
                if(~isfield(optsPhys,'Dmu'))   
                    y0 = [zeros(1,optsPhys.nSpecies);y0];
                end
                return;
          end
            
          %if(~isfield(optsPhys,'HSBulk') && ~isfield(optsNum,'HSBulk'))
          if(~isfield(optsPhys,'Dmu') && ~isfield(optsPhys,'mu'))
              
            % initial guess for mu doesn't really matter
            muInit = zeros(1,optsPhys.nSpecies);

            % however, require relatively good guess for rho

            % rho without interaction ie exp(-(VBack+VAdd )/kBT)
            rhoInit = exp(-(this.Vext+this.VAdd)/optsPhys.kBT);

            % inverse of normalization coefficient
            normalization = repmat( this.IDC.Int*rhoInit./optsPhys.nParticlesS' , size(rhoInit,1) ,1);

            y0 = -optsPhys.kBT*log(normalization) - this.VAdd;
            y0 = [muInit; y0];
          else
             %y0 = -optsPhys.kBT*log(ones(this.IDC.M,1)) - this.Vext;             
           %  rhoInit = exp(-(this.Vext+this.VAdd)/optsPhys.kBT);
            % y0 = - this.VAdd;
             
             y0 = zeros(this.IDC.M,optsPhys.nSpecies);
          end

        end        
        function rho = GetRhoEq(this)
            rho = exp((this.x_eq-this.Vext)/this.optsPhys.kBT);
        end        
        function Preprocess_HardSphereContribution(this)      
            tic
            if(isfield(this.optsNum,'FexNum'))
                fprintf(1,'Computing FMT matrices ...');   
                paramsFex.sigmaS   = this.optsPhys.sigmaS;
                paramsFex.kBT      = this.optsPhys.kBT;            
                paramsFex.Pts      = this.IDC.Pts;
                paramsFex.nSpecies = this.optsPhys.nSpecies;   
                paramsFex.FexNum   = this.optsNum.FexNum;
                
                FexFun             = str2func(['FexMatrices_',this.optsNum.FexNum.Fex]);    
                this.IntMatrFex    = DataStorage(['FexData' filesep class(this.IDC)],FexFun,paramsFex,this.IDC);   
            elseif(isfield(this.optsPhys,'HSBulk') && ~strcmp(this.optsPhys.HSBulk,'Fex_ZeroMap'))
                this.optsNum.FexNum.Fex  = this.optsPhys.HSBulk;                                                               
            else
                this.optsNum.FexNum.Fex  = 'ZeroMap';                           
            end

            fprintf(1,'done.\n');
            t_fex = toc;
            disp(['Fex computation time (sec): ', num2str(t_fex)]);
        end            
        function Preprocess_MeanfieldContribution(this)
            
            if(isfield(this.optsNum,'V2Num'))
                fprintf(1,'Computing mean field convolution matrices ...');   
                paramsFex.V2       = this.optsPhys.V2;
                paramsFex.kBT      = this.optsPhys.kBT;
                paramsFex.FexNum   = this.optsNum.V2Num;
                paramsFex.Pts      = this.IDC.Pts;
                paramsFex.nSpecies = this.optsPhys.nSpecies;   

                FexFun             = str2func(['FexMatrices_',this.optsNum.V2Num.Fex]);    
                this.IntMatrV2     = DataStorage(['FexData' filesep class(this.IDC)],FexFun,paramsFex,this.IDC);   

            elseif(isfield(this.optsPhys,'V2') && ~isfield(this.optsNum,'V2Num'))                
                error('Define V2Num structure in optsNum')
            else
                this.IntMatrV2 = 0;
            end
        end
        function Preprocess_HIContribution(this)
            optsNum  = this.optsNum;
            optsPhys = this.optsPhys;
            
            if(isfield(optsNum,'HINum') && ~isempty(optsNum.HINum))
                this.doHI = true;
            else
                this.doHI = false;
            end

            if(this.doHI)        
                tic
                fprintf(1,'Computing HI matrices ...');   
                paramsHI.optsPhys.HI       = optsPhys.HI;
                paramsHI.optsNum.HINum     = optsNum.HINum;
                paramsHI.optsNum.Pts       = this.IDC.Pts;    
                paramsHI.optsNum.Polar     = 'cart';
                paramsHI.optsPhys.nSpecies = this.optsPhys.nSpecies;
                if(strcmp(optsNum.PhysArea.shape,'HalfSpace_FMT'))
                    this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices_HalfSpace,paramsHI,this.IDC);      
                else
                    this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices2D,paramsHI,this.IDC);      
                end
                fprintf(1,'done.\n');
                t_HI = toc;
                display(['HI computation time (sec): ', num2str(t_HI)]); 
            end
        end
        function Preprocess_ExternalPotential(this)
            PtsCart  = this.IDC.GetCartPts();            
            y1S      = repmat(PtsCart.y1_kv,1,this.optsPhys.nSpecies); 
            y2S      = repmat(PtsCart.y2_kv,1,this.optsPhys.nSpecies);
            ythS     = repmat(this.IDC.Pts.y2_kv,1,this.optsPhys.nSpecies);

            [this.Vext,this.Vext_grad]  = getVBackDVBack(y1S,y2S,this.optsPhys.V1);                  
            if(strcmp(this.IDC.polar,'polar'))
                this.Vext_grad = GetPolarFromCartesianFlux(this.Vext_grad,ythS);                
            end
            this.VAdd  = getVAdd(y1S,y2S,0,this.optsPhys.V1);
        end
        function Preprocess_SubArea(this)
            
            optsNum  = this.optsNum;
            optsPhys = this.optsPhys;

            this.doSubArea = isfield(optsNum,'SubArea');

            if(this.doSubArea)    
                subshapeClass = str2func(optsNum.SubArea.shape);
                this.subArea       = subshapeClass(optsNum.SubArea);
                
                plotSubShape   = optsNum.SubArea;
                plotSubShape.N = [80,80];
                this.subArea.ComputeAll(plotSubShape); 
                
                this.IP            = this.IDC.SubShapePtsCart(this.subArea.GetCartPts());
                this.Int_of_path   = this.subArea.IntFluxThroughDomain(100)*blkdiag(this.IP,this.IP);
            else
                this.Int_of_path   =  zeros(1,2*this.IDC.M);
            end
        end        
        function PlotDynamics(this)
                plotData.optsPhys = this.optsPhys;
                plotData.optsNum  = this.optsNum;
                plotData.data     = this.dynamicsResult;
                
                figure('Position',[0 0 1000 1000]);
                PlotDDFT(plotData);  
        end    
        
        function ResetTemperature(this,T)
            this.optsPhys.kBT = T;
            % Determine saturation point if a 2-Phase system is given
            if(isfield(this.optsPhys,'V2') && ~strcmp(this.optsPhys.HSBulk,'Fex_ZeroMap'))
                [this.optsPhys.rhoGas_sat,...
                     this.optsPhys.rhoLiq_sat,...
                    this.optsPhys.mu_sat,p] = BulkSatValues(this.optsPhys);
                GetCriticalPoint(this.optsPhys);
                % BulkPhaseDiagram(this.optsPhys);
                this.do2Phase = true;
            end
        end
        
        ComputeEquilibrium(this,rho_ig)
        ComputeDynamics(this)
        
        
    end
end