classdef DDFT_2D < Computation
    
    properties (Access = public)                    
        

        IntMatrV2  % integration matrix for mean field two-particle interactions
        IntMatrHI  % integration matrices for hydrodynamic interactions
        IntMatrFex % integration matrix for FMT (hard-sphere) interactions
               
        Vext,Vext_grad,VAdd
        
        x_eq,mu
        FilenameEq,FilenameDyn
        dynamicsResult
        
        doHI,doHIWall,do2Phase
    end
    
    methods (Access = public)          
        function this = DDFT_2D(configuration)      
            if(nargin == 0)
                 configuration = [];
            end             
            this@Computation(configuration);            
            
            if(~isfield(this.optsPhys,'D0') && isfield(this.optsPhys,'mS') && isfield(this.optsPhys,'gammaS'))
                this.optsPhys.D0 = 1./(this.optsPhys.gammaS.*this.optsPhys.mS);
            elseif(~isfield(this.optsPhys,'D0'))
                this.optsPhys.D0 = 1;
            end                
            
            
            if((isfield(this.optsPhys,'Inertial') && this.optsPhys.Inertial) || ...
                     (isfield(this.optsNum,'Inertial') && this.optsNum.Inertial)) 
                 this.optsPhys.Inertial = true;
            else
                this.optsPhys.Inertial = false;
            end
            
        end        
        function res = Preprocess(this)

            Preprocess@Computation(this);
                        
            optsPhys = this.optsPhys;
           
            % Determine hard-sphere contribution to bulk free energy
            if(~isfield(this.optsPhys,'HSBulk'))
                if(isfield(this.optsNum,'FexNum'))
                    this.optsPhys.HSBulk = ['FexBulk_',this.optsNum.FexNum.Fex];
                else
                    this.optsPhys.HSBulk = 'Fex_ZeroMap';
                end
            end                         
            
             % Determine number of species
            if(isfield(optsPhys,'nParticlesS'))
                this.optsPhys.nSpecies = length(optsPhys.nParticlesS); 
            elseif(isfield(optsPhys,'Dmu'))
                this.optsPhys.nSpecies = length(optsPhys.Dmu);
            elseif(isfield(optsPhys,'mu'))
                this.optsPhys.nSpecies = length(optsPhys.mu);
            elseif(isfield(optsPhys,'nSpecies'))
                this.optsPhys.nSpecies = optsPhys.nSpecies;
            end
            
            % Determine saturation point if a 2-Phase system is given
            if(isfield(this.optsPhys,'V2') && ~strcmp(this.optsPhys.HSBulk,'Fex_ZeroMap') && (this.optsPhys.nSpecies == 1))
                [this.optsPhys.rhoGas_sat,...
                     this.optsPhys.rhoLiq_sat,...
                    this.optsPhys.mu_sat,this.optsPhys.p] = BulkSatValues(this.optsPhys);
                GetCriticalPoint(this.optsPhys);
                % BulkPhaseDiagram(this.optsPhys);
                this.do2Phase = true;                                
            end
            if(isfield(this.optsPhys,'mu_sat') && isfield(this.optsPhys,'Dmu'))
                this.mu =  this.optsPhys.mu_sat + this.optsPhys.Dmu;       
            end
            if(isfield(this.optsPhys,'mu'))
                this.mu =  this.optsPhys.mu;
            end            

            res_conv = Preprocess_MeanfieldContribution(this);
            res_fex  = Preprocess_HardSphereContribution(this);
            Preprocess_HIContribution(this);
            Preprocess_ExternalPotential(this);
            Preprocess_SubArea(this);
            
            res = mergeStruct(res_conv,res_fex);
            if(isempty(res))
                res.empty = true;
            end
        end                
        function y0  = getInitialGuess(this,rho_ig)
            
          optsNum  = this.optsNum;
          optsPhys = this.optsPhys;          
          
          if(nargin >= 2)
                if(isscalar(rho_ig))
                    rho_ig = rho_ig*ones(this.IDC.M,this.optsPhys.nSpecies);
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
        function ResetTemperature(this,T)
            if(nargin > 1)
                this.optsPhys.kBT = T;
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
        end        
        function x_ic = GetInitialCondition(this)
            if(~isfield(this.optsPhys,'ModifyEq_to_IC'))
                x_ic = this.x_eq;
            elseif(strcmp(this.optsPhys.ModifyEq_to_IC.Mode,'expY2'))
                %needed elements of struct ModifyEq_to_IC:
                % - Mode {expY2}
                % - a0,a1                
                opts                  = this.optsPhys.ModifyEq_to_IC;
                
                ptsCart_Shifted       = this.IDC.GetCartPts();
                y2C                   = ptsCart_Shifted.y2_kv;
                ptsCart_Shifted.y1_kv = ptsCart_Shifted.y1_kv + opts.a0*exp(-y2C/opts.a1);
                IP                    = this.IDC.SubShapePtsCart(ptsCart_Shifted);
                x_ic                  = IP*this.x_eq;
                
            end
        end        
        
        function res = Preprocess_HardSphereContribution(this)      
            res = struct();
            if(isfield(this.optsNum,'FexNum'))
                fprintf(1,'Computing FMT matrices ...\n');   
                paramsFex.sigmaS   = this.optsPhys.sigmaS;                
                paramsFex.physArea = this.optsNum.PhysArea;                
                paramsFex.nSpecies = this.optsPhys.nSpecies;   
                paramsFex.FexNum   = this.optsNum.FexNum;
                FexFun             = str2func(['FexMatrices_',this.optsNum.FexNum.Fex]);    
                folder             = ['FexData' filesep class(this.IDC) filesep func2str(FexFun)];
                
                paramsFex.Pts      = this.IDC.Pts;
                ClearPts(folder,FexFun,paramsFex,{'physArea','kBT'});                
                paramsFex          = rmfield(paramsFex,'Pts');  
                this.IntMatrFex    = DataStorage(folder,FexFun,paramsFex,this.IDC,[],{'kBT'});   
                
                disp('*** Test FMT matrices ***');
                if(isfield(this.optsNum.FexNum,'Fex') && strcmp(this.optsNum.FexNum.Fex,'FMTRosenfeld_3DFluid'))
                    res = CheckAverageDensities_Rosenfeld_3D(this.IDC,this.IntMatrFex);                    
                end
                
            elseif(isfield(this.optsPhys,'HSBulk') && ~strcmp(this.optsPhys.HSBulk,'Fex_ZeroMap'))
                this.optsNum.FexNum.Fex  = this.optsPhys.HSBulk;                                                               
            else
                this.optsNum.FexNum.Fex  = 'ZeroMap';                           
            end
            fprintf(1,'done.\n');            
        end            
        function res = Preprocess_MeanfieldContribution(this)

            if(isfield(this.optsNum,'V2Num') && ~isempty(this.optsNum.V2Num))
                fprintf(1,'Computing mean field convolution matrices ...\n');   
                if(isfield(this.optsPhys,'r_cutoff'))
                    paramsFex.r_cutoff = this.optsPhys.r_cutoff;
                end
                paramsFex.V2       = this.optsPhys.V2;                                
                paramsFex.FexNum   = this.optsNum.V2Num;
                paramsFex.physArea = this.optsNum.PhysArea;
                paramsFex.Pts      = this.IDC.Pts;
                paramsFex.nSpecies = this.optsPhys.nSpecies;   
                FexFun             = str2func(['FexMatrices_',this.optsNum.V2Num.Fex]);                    
                folder             = ['FexData' filesep class(this.IDC) filesep func2str(FexFun)];
                
                ClearPts(folder,FexFun,paramsFex,{'physArea','kBT'});                
                paramsFex          = rmfield(paramsFex,'Pts');                
                this.IntMatrV2     = DataStorage(folder,FexFun,paramsFex,this.IDC,[],{'kBT'});   %true
                
                res = CheckMeanfieldConvolution(this);

            elseif(isfield(this.optsPhys,'V2') && ~isfield(this.optsNum,'V2Num'))                                
                error('Define V2Num structure in optsNum')
            else
                res = struct();
                this.IntMatrV2 = [];
            end
        end
        function Preprocess_HIContribution(this)
            optsNum  = this.optsNum;
            optsPhys = this.optsPhys;
            
            if(isfield(optsNum,'HINum') && isfield(optsNum.HINum,'HI11'))
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
                    if(isfield(optsNum.HINum,'HIWallFull') && optsNum.HINum.HIWallFull)
                        
                        % remove 1-body wall contribution
                        if(isfield(paramsHI.optsNum.HINum,'Wall'))
                            paramsHI.optsNum.HINum = rmfield(paramsHI.optsNum.HINum,'Wall');
                        end
                        
                        this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices_HalfSpace_Wall,paramsHI,this.IDC);%,true);
                    else
                        this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices_HalfSpace,paramsHI,this.IDC);%,true);
                    end
                else
                    this.IntMatrHI     = DataStorage(['HIData' filesep class(this.IDC)],@HIMatrices2D,paramsHI,this.IDC);      
                end
                fprintf(1,'done.\n');
                t_HI = toc;
                display(['HI computation time (sec): ', num2str(t_HI)]); 
            end
            
            if(isfield(optsNum,'HINum') && isfield(optsNum.HINum,'Wall'))
                this.doHIWall = true;
            else
                this.doHIWall = false;
            end
            
            PtsCart  = this.IDC.GetCartPts();
            if(this.doHIWall)
                if(optsPhys.nSpecies>1)
                    error('HI with wall only implemented for one species');
                end
                tic;
                wallHIfn = str2func(optsNum.HINum.Wall);
                [this.IntMatrHI.DWall,this.IntMatrHI.DDWall] = wallHIfn(PtsCart.y1_kv,PtsCart.y2_kv,optsPhys.HI);
                t_HIWall = toc;
                display(['HI wall computation time (sec): ', num2str(t_HIWall)]);
            else
                this.IntMatrHI.DWall  = ones(size([PtsCart.y1_kv;PtsCart.y2_kv]));
                this.IntMatrHI.DDWall = zeros(size([PtsCart.y1_kv;PtsCart.y2_kv]));
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

        res = ComputeEquilibrium(this,rho_ig,optsIn,miscIn)                
        function ComputeDynamics(this)
            if(this.doHIWall)
                ComputeDynamicsWallHI(this);
            elseif(this.optsPhys.Inertial)
                ComputeDynamicsInertia(this);
            else
                ComputeDynamicsOverdamped(this);
            end
        end
        ComputeDynamicsOverdamped(this)
        ComputeDynamicsInertia(this)
        
        function PlotDensityContours(this,rho,optsPlot)
                       
            optsPhys = this.optsPhys;
            if(isfield(this.optsPhys,'rhoGas_sat'))                
                drho = optsPhys.rhoLiq_sat - optsPhys.rhoGas_sat;
                    
                optsPlot.nContours = optsPhys.rhoGas_sat + 0.1*drho;
                optsPlot.linecolor = 'b';
                optsPlot.linestyle = '--';
                this.IDC.plot(rho,'contour',optsPlot);  hold on;  

                optsPlot.nContours = optsPhys.rhoGas_sat + 0.5*drho;
                optsPlot.linecolor = [0 0.75 0];
                this.IDC.plot(rho,'contour',optsPlot);  hold on;  

                optsPlot.nContours = optsPhys.rhoGas_sat + 0.9*drho;
                optsPlot.linecolor = 'r';
                this.IDC.plot(rho,'contour',optsPlot);  hold on;  
            else
                rho_tMin  = min(min(min(rho_t)));
                rho_tMax  = max(max(max(rho_t)));            
                optsPlot = struct('nContours',rho_tMin + [0.3,0.5,0.7]*(rho_tMax-rho_tMin));
                this.IDC.plot(rho,'contour',optsPlot); 
            end
            
        end        
        function PlotDynamicValue(this,varName,opts)
            if(nargin < 3)
                opts = {};
            end
            
            plotTimes = this.dynamicsResult.t;
            rho_t     = this.dynamicsResult.rho_t;  
            
            if(iscell(varName))
                val       = this.dynamicsResult.(varName{1});
                val2      = this.dynamicsResult.(varName{2});
                name      = [varName{1},'_',varName{2}];
            else
                val       = this.dynamicsResult.(varName);
                val2      = [];
                name      = varName;
            end
                        
            
            figure('Color','white','Position', [0 0 800 800]);
            
            T_n_Max = length(plotTimes);            
            fileNames = [];
            maxV2     = max(max(max(val2)));
            
            for i=1:T_n_Max
            
                t       = plotTimes(i);
                hold off;
                for iSpecies=1:size(val,2)   
                    if(~isempty(val2))                        
                        this.IDC.plot(val2(:,iSpecies,i),'color'); hold on;
                        m = max(max(val2(:,iSpecies,i)));
                        if(m==0)
                            colormap(b2r(0,1));
                        else
                            colormap(b2r(0,m));
                        end
                        colorbar;
                        %set(gca, 'CLim', [min(min(min(val2))) max(max(max(val2)))]);
                    end
                    
                    maxVal = max(max(max(abs(val))));
                    if(size(val,1) == 2*this.IDC.M)                   
                        PlotDensityContours(this,rho_t(:,iSpecies,i));  hold on;                        
                        this.IDC.plotFlux(val(:,iSpecies,i),[],maxVal,1.2,'k'); 
                        %this.IDC.plot(val(:,iSpecies,i),{'flux'});%,[],maxVal,1.5,'k'); 
                    else
                        %this.IDC.plot(val(:,iSpecies,i),'SC');                    
                        this.IDC.plot(val(:,iSpecies,i),'contour');
                    end
                    
                end                
                title(['t = ', num2str(t,'%10.1f')]);               
                set(gca,'fontsize',20);
                set(gca,'linewidth',1.5);
                %xlabel('$y_1$','Interpreter','Latex','fontsize',25);
                %ylabel('$y_2$','Interpreter','Latex','fontsize',25);
                %view([2,5,2]);

                h = get(gca,'xlabel'); set(h,'fontsize',35);
                h = get(gca,'ylabel'); set(h,'fontsize',35);
                h = get(gca,'title');  set(h,'fontsize',35);
                pause(0.2);                
                
                if(IsOption(opts,'save'))
                    fileName = getPDFMovieFile('Movie1',i);
                    disp(fileName);                    
                    save2pdf(fileName,gcf);
                    SaveFigure(fileName(1:end-4));
                    fileNames = [fileNames,' ',fileName];
                end                
                
            end
            
            if(IsOption(opts,'save'))                
                allPdfFiles = [this.FilenameDyn,'_',name,'.pdf'];                
    
                system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
                disp(['Concatenated pdf file saved in: ',allPdfFiles]);            
                system(['del ',fileNames]);           
            end
        end 
        function PlotDynamicValueLine(this,y1P,y2P,varName,opts)
            if(nargin < 5)
                opts = {};
            end
            
            plotTimes = this.dynamicsResult.t;
            rho_t     = this.dynamicsResult.rho_t;  
            
            
            if(~iscell(varName))
                varName = {varName};
            end
            name = [];
            for i = 1:length(varName)
                if(strcmp(varName{i},'U_t'))
                    val{i}    = this.dynamicsResult.UV_t(1:end/2,:,:);
                elseif(strcmp(varName{i},'V_t'))
                    val{i}    = this.dynamicsResult.UV_t(1+end/2:end,:,:);   
                else
                    val{i}    = this.dynamicsResult.(varName{i});                
                end
                name      = [name,'_',varName{1}];
            end                                       
            
            figure('Color','white','Position', [0 0 800 800]);
            
            T_n_Max = length(plotTimes);            
            fileNames = [];            
            
            for i=1:T_n_Max
            
                t       = plotTimes(i);
                hold off;
                

                
                 for iSpecies=1:size(val,2)   
                     
                    for k = 1:length(val)
                        this.IDC.plotLine(y1P,y2P,val{k}(:,iSpecies,i));
                    end
%                     if(~isempty(val2))                        
%                         this.IDC.plot(val2(:,iSpecies,i),'color'); hold on;
%                         colormap(b2r(0,max(max(val2(:,iSpecies,i)))));
%                         colorbar;
%                         %set(gca, 'CLim', [min(min(min(val2))) max(max(max(val2)))]);
%                     end
%                     
%                     maxVal = max(max(max(abs(val))));
%                     if(size(val,1) == 2*this.IDC.M)                   
%                         PlotDensityContours(this,rho_t(:,iSpecies,i));  hold on;                        
%                         this.IDC.plotFlux(val(:,iSpecies,i),[],maxVal,1.2,'k'); 
%                         %this.IDC.plot(val(:,iSpecies,i),{'flux'});%,[],maxVal,1.5,'k'); 
%                     else
%                         %this.IDC.plot(val(:,iSpecies,i),'SC');                    
%                         this.IDC.plot(val(:,iSpecies,i),'contour');
%                     end
%                     
                 end 
                
                title(['t = ', num2str((t))]);               
                set(gca,'fontsize',20);
                set(gca,'linewidth',1.5);
                %xlabel('$y_1$','Interpreter','Latex','fontsize',25);
                %ylabel('$y_2$','Interpreter','Latex','fontsize',25);
                %view([2,5,2]);

                h = get(gca,'xlabel'); set(h,'fontsize',35);
                h = get(gca,'ylabel'); set(h,'fontsize',35);
                h = get(gca,'title');  set(h,'fontsize',35);
                pause(0.2);
                
                
                if(IsOption(opts,'save'))
                    fileName = getPDFMovieFile('Movie1',i);
                    save2pdf(fileName,gcf);
                    fileNames = [fileNames,' ',fileName];
                end                
                
            end
            
            if(IsOption(opts,'save'))                
                allPdfFiles = [this.FilenameDyn,'_',name,'.pdf'];                
    
                system(['C:\pdftk.exe ', fileNames ,' cat output ',allPdfFiles]);    
                disp(['Concatenated pdf file saved in: ',allPdfFiles]);            
                system(['del ',fileNames]);           
            end
        end 
        
        function fig_h = PlotDynamics(this,rec)
            if(nargin == 1)
                rec = false;
            end                                                
            
            plotData.subArea  = this.subArea;
            plotData.optsPhys = this.optsPhys;
            plotData.optsNum  = this.optsNum;
            plotData.data     = this.dynamicsResult;
            plotData.filename = this.FilenameDyn;
            plotData.data.shape = this.IDC;

            figure('Position',[0 0 1000 1000]);
            fig_h = PlotDDFT(plotData,rec);
        end              
        function PostprocessDynamics(this)                    
            optsPhys      = this.optsPhys;
            subArea       = this.subArea;%dynamicsResult.Subspace.subArea;
            rho_t         = this.dynamicsResult.rho_t;
            no_times      = length(this.dynamicsResult.t);
            accFlux       = this.dynamicsResult.Subspace.accFlux;
            nSpecies      = size(rho_t,2);
            rho_ic        = rho_t(:,:,1);
            
            IP            = this.IDC.SubShapePtsCart(subArea.GetCartPts());
            Int_SubOnFull = subArea.ComputeIntegrationVector()*IP;
            
            massError     = zeros(no_times,nSpecies);
                                    
            for iSpecies=1:nSpecies            
                rho          = permute(rho_t(:,iSpecies,:),[1 3 2]);
                rho_diff     = rho-rho_ic(:,iSpecies)*ones(1,no_times);
                massError    = Int_SubOnFull*rho_diff+accFlux(:,iSpecies)';                    
                mass         = (Int_SubOnFull*rho);              
                massErrorRel = massError./mass;  
            end
            
            this.dynamicsResult.Subspace.mass         = mass;
            this.dynamicsResult.Subspace.massError    = massError;            
            this.dynamicsResult.Subspace.massErrorRel = massErrorRel;    
            
            %Add entropy production
            if(isfield(optsPhys,'viscosity') && optsPhys.Inertial && (nSpecies == 1))
            
                entropy = zeros(size(rho_t));
                Diff = this.IDC.Diff;
                
                for nt = 1:size(rho_t,3)
                    %go through time steps
                    
                    rho  = rho_t(:,:,nt);
                    
                    zeta = BulkViscosity(rho,optsPhys);
                    eta  = ShearViscosity(rho,optsPhys);
                    u    = this.dynamicsResult.UV_t(1:end/2,:,nt);
                    v    = this.dynamicsResult.UV_t(1+end/2:end,:,nt);

                    divUV = Diff.div*[u;v];
                    entropy(:,:,nt) = zeta.*(divUV.^2) +...
                                      eta/2.*( ((2*Diff.Dy1)*u - 2/3*divUV).^2 ...
                                             + ((2*Diff.Dy2)*v - 2/3*divUV).^2 ...
                                             + 2*(Diff.Dy1*v+Diff.Dy2*u).^2 );     
                end                
                this.dynamicsResult.entropy = entropy;
            end
            
            
        end
    end
end