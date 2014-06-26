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
        
        doHI,doSubArea
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
            
            shape        = optsNum.PhysArea;
            if(isfield(optsNum,'V2Num'))
                shape.Conv = optsNum.V2Num;
            else                
                shape.Conv = 0;
            end
            
            shapeClass = str2func(optsNum.PhysArea.shape);
            this.IDC   = shapeClass(shape);            

            this.IDC.ComputeAll(optsNum.PlotArea); 
            PtsCart                = this.IDC.GetCartPts();

            if(isfield(optsPhys,'nParticlesS'))
                this.optsPhys.nSpecies = length(optsPhys.nParticlesS); 
            elseif(isfield(optsPhys,'Dmu'))
                this.optsPhys.nSpecies = length(optsPhys.Dmu);
                [this.optsPhys.rhoGas_sat,...
                 this.optsPhys.rhoLiq_sat,...
                 this.optsPhys.mu_sat,p] = BulkSatValues(this.optsPhys);
            end
            
            if(~isfield(optsPhys,'sigmaS'))                
                this.optsPhys.sigmaS = ones(this.optsPhys.nSpecies,1);
            end
                        
            %************************************************
            %****************  Preprocess  ****************
            %************************************************  

            tic
            
            % Hard sphere contribution             
            if(isfield(optsNum,'FexNum'))
                fprintf(1,'Computing FMT matrices ...');   
                paramsFex.V2       = optsPhys.V2;
                paramsFex.kBT      = optsPhys.kBT;            
                paramsFex.Pts      = this.IDC.Pts;
                paramsFex.nSpecies = this.optsPhys.nSpecies;   
                paramsFex.FexNum   = optsNum.FexNum;
                
                FexFun             = str2func(['FexMatrices_',optsNum.FexNum.Fex]);    
                this.IntMatrFex    = DataStorage(['FexData' filesep class(this.IDC)],FexFun,paramsFex,this.IDC);   
            elseif(isfield(optsPhys,'HSBulk'))
                this.optsNum.FexNum.Fex  = optsPhys.HSBulk;                                                               
            else
                this.optsNum.FexNum.Fex  = 'ZeroMap';                           
            end

            fprintf(1,'done.\n');
            t_fex = toc;
            disp(['Fex computation time (sec): ', num2str(t_fex)]);
                        
            if(isfield(optsNum,'V2Num'))
                fprintf(1,'Computing mean field convolution matrices ...');   
                paramsFex.V2       = optsPhys.V2;
                paramsFex.kBT      = optsPhys.kBT;
                paramsFex.FexNum   = optsNum.V2Num;
                paramsFex.Pts      = this.IDC.Pts;
                paramsFex.nSpecies = this.optsPhys.nSpecies;   

                FexFun             = str2func(['FexMatrices_',optsNum.V2Num.Fex]);    
                this.IntMatrV2     = DataStorage(['FexData' filesep class(this.IDC)],FexFun,paramsFex,this.IDC);   

                fprintf(1,'done.\n');
                t_fex = toc;
                disp(['Fex computation time (sec): ', num2str(t_fex)]);
            elseif(isfield(this.optsPhys,'V2') && ~isfield(this.optsNum,'V2Num'))                
                error('Define V2Num structure in optsNum')
            else
                this.IntMatrV2 = 0;
            end

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
                this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices2D,paramsHI,this.IDC);      
                fprintf(1,'done.\n');
                t_HI = toc;
                display(['HI computation time (sec): ', num2str(t_HI)]); 
            end
            
            % External Potential
            y1S     = repmat(PtsCart.y1_kv,1,this.optsPhys.nSpecies); 
            y2S     = repmat(PtsCart.y2_kv,1,this.optsPhys.nSpecies);
            ythS     = repmat(this.IDC.Pts.y2_kv,1,this.optsPhys.nSpecies);

            [this.Vext,this.Vext_grad]  = getVBackDVBack(y1S,y2S,this.optsPhys.V1);                  
            if(strcmp(this.IDC.polar,'polar'))
                this.Vext_grad = GetCartesianFromPolarFlux(this.Vext_grad,ythS);
            end
            this.VAdd  = getVAdd(y1S,y2S,0,this.optsPhys.V1);

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
        function y0 = getInitialGuess(this)
            
          optsNum  = this.optsNum;
          optsPhys = this.optsPhys;          
            
          %if(~isfield(optsPhys,'HSBulk') && ~isfield(optsNum,'HSBulk'))
          if(~isfield(optsPhys,'Dmu'))
              
            % initial guess for mu doesn't really matter
            muInit=zeros(1,optsPhys.nSpecies);

            % however, require relatively good guess for rho

            % rho without interaction ie exp(-(VBack+VAdd )/kBT)
            rhoInit = exp(-(this.Vext+this.VAdd)/optsPhys.kBT);

            % inverse of normalization coefficient
            normalization = repmat( this.IDC.Int*rhoInit./optsPhys.nParticlesS' , size(rhoInit,1) ,1);

            y0 = -optsPhys.kBT*log(normalization) - this.VAdd;
            y0 = [muInit; y0];
          else
             %y0 = -optsPhys.kBT*log(ones(this.IDC.M,1)) - this.Vext;             
             y0 = zeros(this.IDC.M,optsPhys.nSpecies);
          end

        end        
        function rho = GetRhoEq(this)
            rho = exp((this.x_eq-this.Vext)/this.optsPhys.kBT);
        end
        
        ComputeEquilibrium(this,rho_ig)
        ComputeDynamics(this)
        
        
    end
end