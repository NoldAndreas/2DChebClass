classdef ContactLineHS < DDFT_2D
    
    properties (Access = public)                                           
        %(1) 1D Results
        rho1D_lg,rho1D_wl,rho1D_wg,alpha_YCA  
        ST_1D  %1D surface tensions
        
        y1_SpectralLine        
        hII = [], hIII = [],hIV = [], hContour = []
        disjoiningPressure_II=[],disjoiningPressure_IV=[]
        y1_I = [], hI = []
        
        y2=[],Int_y2=[],DiffY2=[]
        
        IsolineInterfaceY2=[]        
        FittedInterfaceY2=[]
               
        %structure with information of the adsorption isotherm        
        AdsorptionIsotherm 
        
        grandPot,grandPot2
        disjoiningPressureCheck
        
        wIntLoc,wIntHS,f_loc,f_hs
        CA_deg_measured,CA_deg_measured_error
        
        saveFigs = true;
    end    
    methods (Access = public)          
        function this = ContactLineHS(configuration)        
             if(nargin == 0)
                 configuration = [];
             else
                configuration.optsNum.PhysArea.shape = 'HalfSpace_FMT';                
             end
             
             configuration.optsNum.PhysArea.L2_AD  = configuration.optsNum.PhysArea.L2;
             configuration.optsNum.PhysArea.y2wall = 0;
             configuration.optsNum.PhysArea.h      = 1;             
             
             if(~isfield(configuration.optsNum.PhysArea,'N2bound'))
                 configuration.optsNum.PhysArea.N2bound = max(10,2*round(configuration.optsNum.PhysArea.N(2)/6));
             end
             
             this@DDFT_2D(configuration);
        end                
        
        %Preprocessing and testing of results
        function res = Preprocess(this)
            global dirDataOrg
            ChangeDirData([dirDataOrg filesep 'deg',num2str(this.optsNum.PhysArea.alpha_deg,3)]);
            res = Preprocess@DDFT_2D(this);     
            TestPreprocess(this);        
            
            if(isfield(this.optsPhys,'rhoLiq_sat') && ...
                    isfield(this.optsPhys,'rhoGas_sat'))
                ComputeST(this);
            end
        end        
        function TestPreprocess(this)                          
                        
            R = this.IDC.R;                                                  
                        
            %(3) Test Integration
            disp('*** Test integration vector ***');
            a = 2;
            
            %(a) Integration in HalfSpace
            if(this.IDC.N1 > 1)
                pts   = this.IDC.GetCartPts();
                r     = sqrt(pts.y1_kv.^2+(pts.y2_kv-R).^2);
                testF = exp(-(r/a).^2);    
                PrintErrorPos(this.IDC.Int*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);    
            end            
            pts   = this.IDC.GetCartPts(zeros(this.IDC.N2,1),this.IDC.Pts.y2);
            [h1,h2,Int_1D] = this.IDC.ComputeIntegrationVector();
            testF = exp(-pts.y2_kv);              
            PrintErrorPos(Int_1D*testF-exp(-min(pts.y2_kv)),'Integration of exp(-y_2)');

            %(b) Integration in Composed Half Space
            pts   = this.IDC.AD.GetCartPts();
            r     = sqrt(pts.y1_kv.^2+pts.y2_kv.^2);
            testF = exp(-(r/a).^2);
            Inth  = this.IDC.AD.ComputeIntegrationVector();
            PrintErrorPos(Inth*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);
            disp('*****');
                                    
        end         
                
        %1D surface tensions
        function [om,rho1D,params] = Compute1D(this,WLWGLG)
                        
            optss            = this.optsPhys;   
            Fex_Num          = this.optsNum.FexNum;                        
            
            if(isnumeric(WLWGLG))
                optss.eta        = WLWGLG;
                optss.rho_iguess = WLWGLG*6/pi;
                [rho1D,params]   = FMT_1D(this.IDC,this.IntMatrFex,optss,Fex_Num,[]);   
            elseif(strcmp(WLWGLG,'WL'))
                optss.rho_iguess = this.optsPhys.rhoLiq_sat;
                [rho1D,params] = FMT_1D(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv,{}); %'plot'
                this.ST_1D.om_wallLiq = params.Fex;
                this.rho1D_wl         = rho1D;                
            elseif(strcmp(WLWGLG,'WG'))
                optss.rho_iguess = this.optsPhys.rhoGas_sat;
                [rho1D,params] = FMT_1D(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv,{});%'plot'
                
                this.ST_1D.om_wallGas = params.Fex;
                this.rho1D_wg         = rho1D;                
            elseif(strcmp(WLWGLG,'LG') || strcmp(WLWGLG,'LV'))
                rhoLiq_sat       = this.optsPhys.rhoLiq_sat;
                rhoGas_sat       = this.optsPhys.rhoGas_sat;
                Pts              = this.IDC.Pts;
                theta_CS         = this.IDC.alpha;  
            
                optss.rho_iguess = (rhoLiq_sat+rhoGas_sat)/2 + ...
                                  (rhoLiq_sat-rhoGas_sat)/2*tanh(Pts.y1*sin(theta_CS));
                [rho1D,params] = FMT_1D_Interface(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv,false);                
                
                this.ST_1D.om_LiqGas  = params.Fex;
                this.rho1D_lg         = rho1D;
            end
            
            om  = params.Fex;                                  
            fprintf(['Surface tension ',WLWGLG,' = ',num2str(om),'\n']);
        end                
        function alpha_YCA = ComputeST(this)                        
            Compute1D(this,'WL');            
            Compute1D(this,'WG');
            if(this.IDC.N1 > 1)
                Compute1D(this,'LG');                      
                alpha_YCA = ComputeContactAngle(this.ST_1D.om_wallGas,...
                            this.ST_1D.om_wallLiq,...
                            this.ST_1D.om_LiqGas);            
                this.alpha_YCA = alpha_YCA;
            else
                alpha_YCA = [];
            end            
        end        
        
        % Initialization
        function InitInterpolation(this,bounds1,bounds2)
            PlotAreaCart = struct('y1Min',bounds1(1),'y1Max',bounds1(2),...
                      'y2Min',bounds2(1),'y2Max',bounds2(2),...
                      'N1',100,'N2',100);
            this.IDC.InterpolationPlotCart(PlotAreaCart,true);                  
        end              
        function SetV1(this,epw)
            this.optsPhys.V1.epsilon_w = epw;
            PtsCart                    = this.IDC.GetCartPts();
            this.VAdd                  = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,this.optsPhys.V1);
        end       
        function SetkBT(this,kBT)
            this.optsPhys.kBT = kBT;
            this.Preprocess_BulkValues(); 
        end
        
        function y0  = getInitialGuess(this,rho_ig)
                    
        	rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
            rhoGas_sat    = this.optsPhys.rhoGas_sat;
            kBT           = this.optsPhys.kBT;            
            
            if(nargin < 2)
                if(~isempty(this.rho1D_lg))
                    p = (this.rho1D_lg-rhoGas_sat)/(rhoLiq_sat-rhoGas_sat);    
                else
                    p = 1;
                end
                rho_ig    = kron(p,this.rho1D_wl) + kron(1-p,this.rho1D_wg);         
            else
                if(length(rho_ig)==1)
                    rho_ig = rho_ig * ones(this.IDC.M,1);                
                end
            end
            y0        = kBT*log(rho_ig)+this.Vext;
        end
        
        %Plot functions        
        PlotContourResults(this,plain)
        function PlotDensityResult(this,DP)                                    
            rho      = GetRhoEq(this);
            PlotArea = this.optsNum.PlotAreaCart;
            
            figure('Color','white','Position',[0 0 800 800]);
            this.IDC.plot(rho);
            zlabel('$n\sigma^3$','Interpreter','Latex','fontsize',26);
            if(isfield(PlotArea,'zMax'))
                zlim([0 PlotArea.zMax]);
            end
            colormap(hsv);
            set(gca, 'CLim', [0, 1.0]);            
            pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 5]);
            view([-10 5 3]);
            
            SaveCurrentFigure(this,'Density');   
            if((nargin == 1) || ~strcmp(DP,'DP'))                
                return;
            end
            
            %****** Plot disjoining pressures
            %f2         = figure;
            y1         = this.y1_SpectralLine.Pts.y;
            My1        = length(y1);
            baseline_z = 7;
            factor_z   = 40;            
            hold on;           
            YL = 0.5*ones(My1,1);
            
            plot3(y1,YL,baseline_z*ones(My1,1),'k','linewidth',1.5);
            
           % if((this.optsNum.PhysArea.alpha_deg<90) || (this.optsNum.PhysArea.alpha_deg>=120))
            %if(this.optsNum.PhysArea.alpha_deg ~= 90)
                [~,DeltaY1_III] = this.ComputeDeltaFit();
          %  else
           %     DeltaY1_III = 0;
           % end
                plot3(this.y1_I+DeltaY1_III,0.5*ones(size(this.y1_I)),...
                        baseline_z+factor_z*GetDisjoiningPressure_I_ell(this,this.hI),...
                        'b','linewidth',1.5);
           % end
            
            plot3(y1,YL,baseline_z+factor_z*this.disjoiningPressure_II,...
                'r','linewidth',1.5);     
            
            if(this.optsNum.PhysArea.alpha_deg ~= 90)
                markIII = (this.hIII < PlotArea.y2Max);
                Pi_III = GetDisjoiningPressure_III(this);            
                plot3(y1(markIII),YL(markIII),baseline_z+factor_z*Pi_III(markIII),'k','linewidth',1.5);
            end
            
            plot3(y1,YL,baseline_z+factor_z*this.disjoiningPressure_IV,...
                'g','linewidth',1.5);     
            
             zlim([0 baseline_z+1]);  
             set(gca,'ZTick',[0,baseline_z]);
             
             %labelDP = '$\Pi \sigma^3/\varepsilon$';
             %plot::Text3d(labelDP,[max(y1),0,baseline_z]);%,'Interpreter','Latex','fontsize',25);
            
    

            SaveCurrentFigure(this,'3D_DP');                        
        end
        function PlotInterfaceResults(this)
            %*******************************************************
            % ***************** Interface Plots ********************
            %*******************************************************
            figure('Color','white','Position',[0 0 800 1000],'name','1D Interface Plots');
            subplot(3,1,1);
            this.IDC.do1DPlotParallel(this.rho1D_lg); 
            title('Liquid-Gas Interface');
            ylabel('$\varrho$','Interpreter','Latex');
            xlabel('$y_1$','Interpreter','Latex');

            subplot(3,1,2);
            this.IDC.do1DPlotNormal(this.rho1D_wg);
            title('Wall-Gas Interface');
            ylabel('$\varrho$','Interpreter','Latex');  
            xlabel('$y_2$','Interpreter','Latex');

            subplot(3,1,3);
            this.IDC.do1DPlotNormal(this.rho1D_wl);
            title('Wall-Liquid Interface');
            ylabel('$\varrho$','Interpreter','Latex'); 
            xlabel('$y_2$','Interpreter','Latex');
            
            SaveCurrentFigure(this,'Interfaces');
        end
        function f2 = PlotDisjoiningPressures(this,opts)
            
            if(nargin == 1)
                opts = {};
            end

            y1   = this.y1_SpectralLine.Pts.y;    

            f2 = figure('Color','white','Position',[0 0 1000 750]);            
            plot([-10 35],[0 0],'k'); hold on;
            if(~isempty(this.disjoiningPressure_II))
                plot(y1,this.disjoiningPressure_II,'k--');
            end
            if(~isempty(this.disjoiningPressure_IV) && ~ IsOption(opts,'noDP_IV'))
                plot(y1,this.disjoiningPressure_IV,'g--');
            end    
            
            if( abs(90 - this.alpha_YCA*180/pi) > 5)
                markFT = (this.hIII < 25);
                hhFT   = GetDisjoiningPressure_III(this);
                plot(y1(markFT),hhFT(markFT),'k');
            end

            [DeltaY1_II,DeltaY1_III] = this.ComputeDeltaFit();
            dP1D                     = GetDisjoiningPressure_I_ell(this,this.hI);
            %plot(this.y1_I+DeltaY1_II,dP1D,'k-.','linewidth',1.5);
            plot(this.y1_I+DeltaY1_III,dP1D,'k-.');


            xlim([min(y1) max(y1)]);
            %ylim([1.1*min(dp1D) 3*max(0.001,max(dp1D))]);%[-0.035 0.005])
            set(gca,'fontsize',20);
            xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
            ylabel('$\Pi \sigma^3/\varepsilon$','Interpreter','Latex','fontsize',25);
            %set(gca,'YTick',-0.1:0.01:0);    
             
            SaveCurrentFigure(this,'DisjoiningPressures');
            %SaveFigure(this,'DisjoiningPressures');

        end
        PlotDensitySlicesNormalInterface(this,y1Lim)
        PlotDensitySlices(this)
        PlotDensitySlicesMovie(this);  
        PlotInterfaceFittingQuality(this,nts);
                
        %Adsorption Isotherm
        ComputeAdsorptionIsotherm(this,n,drying)
        FittingAdsorptionIsotherm(this,FT_Int,n)
        [theta_sumrule,errTheta,I,errI] = SumRule_AdsorptionIsotherm(this,ST_LG)
        function [Pi_I,V_I] = GetDisjoiningPressure_I(this)
            rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
            rhoGas_sat     = this.optsPhys.rhoGas_sat;    
            mu_sat         = this.optsPhys.mu_sat;
            
            if(IsDrying(this)) %drying
                Pi_I = (this.AdsorptionIsotherm.mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);
            else %wetting
                Pi_I = -(this.AdsorptionIsotherm.mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);
            end
            
            fT   = abs(this.AdsorptionIsotherm.FT);
            IntM = zeros(length(fT));
            vh   = zeros(1,length(fT));
            for i = 2:length(fT)        
                %Matrix for integral int(v(y),y=y1..Y)
                h           = fT(i) - fT(i-1);
                vh([i-1,i]) = vh([i-1,i]) + h/2;
                IntM(i,:)   = vh;
            end                            
            V_I = -IntM*Pi_I;                
        end        

        function dP1D = GetDisjoiningPressure_I_ell(this,ell) 
            fT   = abs(this.AdsorptionIsotherm.FT);
            Pi_I = GetDisjoiningPressure_I(this);
            
            ellAD = min(fT)+ell;
            IP    = barychebevalMatrix(fT,ellAD);
            dP1D  = IP*Pi_I;
            dP1D(ellAD<min(fT)) = 0;
            dP1D(ellAD>max(fT)) = 0;
        end        
        function [Pi_III] = GetDisjoiningPressure_III(this)
            D       = this.y1_SpectralLine.Diff.Dy;
            D2      = this.y1_SpectralLine.Diff.DDy;
                 
            Pi_III  = -this.ST_1D.om_LiqGas*(D2*this.hIII)./((1+(D*this.hIII).^2).^1.5);            
        end
        function drying = IsDrying(this)
            %if(~isempty(this.AdsorptionIsotherm))
%                drying = (min(this.AdsorptionIsotherm.FT) < 0);
%            else
                drying = (this.optsNum.PhysArea.alpha_deg > 90);
%            end
        end
        
        %Compute functions 
        function res = ComputeEquilibrium(this,optsIn)
            %if(isfield(this.optsNum,'maxComp_y2'))
                optsIn.maxComp_y2 = this.optsNum.maxComp_y2;
                miscIn.mark       = (this.IDC.GetCartPts.y2_kv <= ...
                                        this.optsNum.maxComp_y2);
                res = ComputeEquilibrium@DDFT_2D(this,[],optsIn,miscIn);
            %else
%                ComputeEquilibrium@DDFT_2D(this,this.optsPhys.rhoLiq_sat,optsIn,[]);
%            end
            
        end
        sol = Compute(this)                     
        
        %Compute disjoining pressure
        Compute_DisjoiningPressure_II(this,y1Int)
        Compute_DisjoiningPressure_IV(this)
        [errRel,estTheta,error_estTheta,I,err] = SumRule_DisjoiningPressure(this,II_or_IV)      
        [err,eps,pi_II] = SumRuleIIError(this,interval_y1)
        function d = Compute_LiquidVapourInterfaceThickness(this)
            
            %Find y1 where (n - n_vap)/(n_liq-n_vap) = 0.05/0.95
            nLV = this.rho1D_lg;
            x1  = this.IDC.Pts.x1;   
            nV  = this.optsPhys.rhoGas_sat;
            nL  = this.optsPhys.rhoLiq_sat;
            
            pTarget = 0.95; xP = fsolve(@GetN_LV,0);
            pTarget = 0.05; xM = fsolve(@GetN_LV,0);
            
            d = (this.IDC.PhysSpace1(xP) - this.IDC.PhysSpace1(xM))*sin(this.IDC.alpha);
            
            function err = GetN_LV(x)
                IP  = barychebevalMatrix(x1,x);
                err = (IP*nLV - nV)/(nL-nV)-pTarget;
            end
        end
        
        %Compute height profiles
        Compute_hI(this);
        Compute_hII(this,II_or_IV);
        hIII = Compute_hIII(this,rho);
        function Compute_hContour(this,level)
            rho_eq       = this.GetRhoEq(); 
            N            = this.y1_SpectralLine.N;
            drho         = (this.optsPhys.rhoLiq_sat - this.optsPhys.rhoGas_sat);

            if((nargin < 2) || isempty(level))
                level = 0.5;
            end
            rhoV         = this.optsPhys.rhoGas_sat + level*drho;

            fsolveOpts   = optimset('Display','off');
            y2Cart       = zeros(N,1);

            y2I = 10;
            for i = 1:N
                pt.y1_kv  = this.y1_SpectralLine.Pts.y(i);        
                y2Cart(i) = fsolve(@rhoX2,y2I,fsolveOpts);                    
                %y2I       = max(y2Cart(i),4);        
                y2I       = y2Cart(i);
            end

            this.hContour = y2Cart;

            function z = rhoX2(y2)
                pt.y2_kv = y2;
                IP = this.IDC.SubShapePtsCart(pt);
                z  = IP*rho_eq-rhoV;
            end    

        end
        function [DeltaY1_II,DeltaY1_III] = ComputeDeltaFit(this)
            
            if(this.optsNum.PhysArea.alpha_deg == 90)
                DeltaY1_II  = 0;
                DeltaY1_III = 0;
                return;
            end

            dP1D            = GetDisjoiningPressure_I_ell(this,this.hI);
            [min_I,i_I]     = min(dP1D);
            [min_III,i_III] = min(GetDisjoiningPressure_III(this));
            DeltaY1_III     = this.y1_SpectralLine.Pts.y(i_III) - this.y1_I(i_I);
            
            hS           = max(this.hI);
            h0           = min(this.hIII);

            fsolveOpts   = optimset('Display','off');            
            f            = this.hIII-h0;
            %f = GetDisjoiningPressure_III(this);
            [~,j]        = max(this.hIII);
         %   [DeltaY1_III,~,exitflag]  = fsolve(@fX,this.y1_SpectralLine.Pts.y(j),fsolveOpts);            
         %   if(exitflag < 1)
         %       cprintf('*r','ComputeDeltaFit: Fitting hIII vs hI: no solution found');
         %   end
            
            f            = this.hII;
            [~,j]        = max(this.hII);
            [DeltaY1_II,~,exitflag] = fsolve(@fX,this.y1_SpectralLine.Pts.y(j),fsolveOpts);            
            if(exitflag < 1)
                cprintf('*r','ComputeDeltaFit: Fitting hII vs hI: no solution found');
            end
            
            function z = fX(y1)                
                IP = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1);
                z  = IP*f-hS;
            end    

        end

        %Postprocess
        [CA_measured,err] = MeasureContactAngle(this,type,yInt)
        [y1Cart]          = ComputeInterfaceContourY2(this,level,y2)
        function fullName = SaveCurrentFigure(this,filename,foldername)
            if(nargin < 3)
                foldername = 'Equilibrium';
            end
            
            if(strcmp(foldername,'Equilibrium') && ~isempty(this.FilenameEq))
                [~,fn]   = fileparts(this.FilenameEq);
                if(~isempty(filename))
                    filename = [fn,'_',filename];            
                else
                    filename = fn;
                end
            end
            if(nargin>2)
                filename = [foldername filesep filename];
            end
            fullName = SaveCurrentFigure@Computation(this,filename);            
        end
        function PostprocessDynamics(this,y2Lim)

            if(nargin < 2)
                y2Lim = [5 7.5];
            end
            
            PostprocessDynamics@DDFT_2D(this);
            
            if(~isfield(this.optsNum,'PlotAreaCart'))
                return;
            end
            
            this.dynamicsResult.pathlines = {};   
            this.dynamicsResult.streaklines = {};   
            for y1_i = -3:3:3
                for y2_i = 0.5:2:6.5
                    this.dynamicsResult.pathlines{end+1}   = this.GetPathlines(y1_i,y2_i);
                    this.dynamicsResult.streaklines{end+1} = this.GetStreaklines(y1_i,y2_i);                    
                end
            end
                        
            %Compute Position Of Contact Line
            rho_t       = this.dynamicsResult.rho_t;
            nSpecies    = size(rho_t,2);
            times       = this.dynamicsResult.t;
            no_t        = size(rho_t,3);
                        
            if((nSpecies > 1))
                return;
            end            
            rho0        = 0.5*(this.optsPhys.rhoGas_sat + this.optsPhys.rhoLiq_sat);
                
            pt.y2_kv    = min(this.IDC.GetCartPts.y2_kv);
            y1_0_Cart   = zeros(no_t,1);
            vel_0_Cart  = zeros(no_t,1);
            ca_deg      = zeros(no_t,1);
            y2_Interface= (this.optsNum.PlotAreaCart.y2Min:0.1:this.optsNum.PlotAreaCart.y2Max)';
            %fsolveOpts  = optimset('Display','off');
            
            this.dynamicsResult.fittedInterface.y1 = zeros(length(y2_Interface),nSpecies,no_t);
            this.dynamicsResult.fittedInterface.y2 = zeros(length(y2_Interface),nSpecies,no_t);
            this.dynamicsResult.ak                 = zeros(2,no_t);
            
            ak = [];            
            for nt = 1:no_t
                rho_eq        = rho_t(:,1,nt);
                %y1_0_Cart(nt) = fsolve(@rhoX1,0,fsolveOpts);
                
                [y1_0_Cart(nt),ca_deg(nt),y1_Interface{nt},ak] = ExtrapolateInterface(this,rho_eq,y2_Interface,ak,y2Lim);
                this.dynamicsResult.fittedInterface.y1(:,1,nt) = y1_Interface{nt};
                this.dynamicsResult.fittedInterface.y2(:,1,nt) = y2_Interface;
                this.dynamicsResult.ak(:,nt)                   = ak;
            end         
            
            vel_0_Cart(1) = (y1_0_Cart(2)-y1_0_Cart(1))/(times(2)-times(1));
            for nt = 2:(size(rho_t,3)-1)      
                vel_0_Cart(nt) = (y1_0_Cart(nt+1)-y1_0_Cart(nt-1))/(times(nt+1)-times(nt-1));
            end
            vel_0_Cart(end) = (y1_0_Cart(end)-y1_0_Cart(end-1))/(times(end)-times(end-1));

            this.dynamicsResult.contactlinePos_y1_0 = (y1_0_Cart-y1_0_Cart(1));
            this.dynamicsResult.contactangle_0.val = ca_deg;
            this.dynamicsResult.contactangle_0.str = '\theta';
            this.dynamicsResult.contactlineVel_y1_0 = vel_0_Cart;
            
            if(~isempty(this.dynamicsResult.pathlines{1}))
                this.dynamicsResult.pathlines_MovingFrameOfReference   = this.dynamicsResult.pathlines;
                this.dynamicsResult.streaklines_MovingFrameOfReference = this.dynamicsResult.streaklines;

                for jj = 1:length(this.dynamicsResult.pathlines)                
                    this.dynamicsResult.pathlines_MovingFrameOfReference{jj}.y1 = ...
                        this.dynamicsResult.pathlines{jj}.y1 - this.dynamicsResult.contactlinePos_y1_0;

                    this.dynamicsResult.streaklines_MovingFrameOfReference{jj}.y1 = ...
                        this.dynamicsResult.streaklines{jj}.y1 - this.dynamicsResult.contactlinePos_y1_0;
                end
            end
            
            function z = rhoX1(y1)
                pt.y1_kv = y1;
                IP       = this.IDC.SubShapePtsCart(pt);
                z        = IP*rho_eq-rho0;
            end 
        end
        [contactlinePos,contactAngle_deg,y2Interface,ak] = ExtrapolateInterface(this,rho,y2,ak,y2Lim,opts,t)
        function sliplength = GetSliplength(this,delta_y1,y2Bar)
            pt.y1_kv   = this.dynamicsResult.contactlinePos_y1_0(end) + delta_y1;
            pt.y2_kv   = y2Bar + 0.5;
            
            IP         = this.IDC.SubShapePtsCart(pt);
            
            uBar       = IP*this.dynamicsResult.UV_t(1:end/2,:,end);
            duBar      = IP*this.IDC.Diff.Dy2*this.dynamicsResult.UV_t(1:end/2,:,end);
            sliplength = uBar/duBar - y2Bar;
        end
        
        function PlotSlipEstimate(this)
            
            u     = this.dynamicsResult.UV_t(1:end/2,:,end);
            CLpos = this.dynamicsResult.contactlinePos_y1_0(end);
            
            opts  = struct('plain',true,'CART',true,'dist0',true);
                       
            PlotSub(-3);%In the vapour phase
            PlotSub(0); %At the contact line
            PlotSub(3); %In the liquid phase                        
            
            function PlotSub(deltaY1)
                figure('Position',[0 0 200 200]);
                this.IDC.plotLine(CLpos*[1 1]+deltaY1,[0.5 5.5],u,opts);
                slips = [GetSliplength(this,deltaY1,1),...
                         GetSliplength(this,deltaY1,1.5),...
                         GetSliplength(this,deltaY1,2),...
                         GetSliplength(this,deltaY1,2.5)];
                
                pbaspect([1 1 1]);
                xlabel('$y_2$'); ylabel('$u$');    
                text(2,0,['$\Sliplength = [',num2str(min(slips),3),','...
                                            num2str(max(slips),3),']$']);
                                        
                filename = [this.FilenameDyn,'_','SlipEstimate_DeltaY_',num2str(deltaY1)];
                SaveCurrentFigure(this,filename);                                        
            end
                        
        end
            
        %to delete
        [f,y1]  = PostProcess(this,y1Int)
        I       = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);
        
        PlotDisjoiningPressureAnalysis(this)    
        PlotInterfaceAnalysisY1(this)
        [y2,theta] = PlotInterfaceAnalysisY2(this,yInt)                

        function disjoingPressure1DCheck(this)
            rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
            rhoGas_sat     = this.optsPhys.rhoGas_sat;
            kBT            = this.optsPhys.kBT;
            mu_sat         = this.optsPhys.mu_sat;    
            Dmu            = this.optsPhys.Dmu;    
            R              = this.optsPhys.sigmaS/2; 
            getFex         = str2func(['Fex_',this.optsNum.FexNum.Fex]);    

            optsPhys = this.optsPhys;        

            rho_Bulk       = rhoGas_sat;
            mu             = mu_sat+Dmu;

            rho            = this.rho_eq;    

            %******************************************

            %Get Bulk value and compare:        
            [h1s,f_hs_Bulk,h2s,h3s] = FexBulk_FMTRosenfeld_3DFluid(rho_Bulk,kBT);

            Phi_r             = str2func(optsPhys.V2.V2DV2);        
            [h1s,h2s,alpha]   = Phi_r(0);
            f_attr_Bulk       = alpha*rho_Bulk^2;
            f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
            f_Vmu_Bulk        = -mu*rho_Bulk;

            floc_Bulk         = f_id_Bulk + f_attr_Bulk + f_Vmu_Bulk;

            %Compute excess grand potential        
            f_id           = kBT*rho.*(log(rho)-1);
            [h1s,h2s,f_hs] = getFex(rho,this.IntMatrFex,kBT,R);
            f_attr         = 0.5*rho.*(this.Conv*rho);

            f_Vmu          = rho.*(this.VAdd- mu);        

            f_loc  = f_id + f_attr + f_Vmu;        

            this.f_loc = f_loc;
            this.f_hs  = f_hs;

            %Initialize Integration Boxes which allow for integration normal to
            %wall
            y2Max = this.optsNum.maxComp_y2 + 10;

            %y1 = this.y1;
        %    f             = zeros(size(y1));
        %    f2            = zeros(size(y1));

          %  wIntLoc       = zeros(length(y1),length(rho));
         %   wIntHS        = zeros(length(y1),length(f_hs));

            PtsCart       = this.HS.GetCartPts();   
            [h1,dVAdd]    = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,this.optsPhys.V1);


            [I,w,weights,IP,pts] = this.HS.doIntFLine([inf inf],[0.5 y2Max],[],'CHEB');
            [h1,dVadd_i]         = getVAdd(pts.y1_kv,pts.y2_kv,0,this.optsPhys.V1);
            IP                   = IP(:,this.HS.Pts.y1_kv == inf);
            fB_AI                = zeros(size(this.AdsorptionIsotherm_FT));
            for i = 1:length(this.AdsorptionIsotherm_FT)
                optss              = this.optsPhys;
                optss.Dmu          = this.AdsorptionIsotherm_Mu(i) - mu_sat;
                optss.rho_iguess   = 'WL';
                [rho_wl,postParms] = FMT_1D(this.HS,this.IntMatrFex,optss,this.optsNum.FexNum,this.Conv,false);

                rho_i    = this.AdsorptionIsotherm_rho(i,:);
                fB_AI(i) = - weights*(dVadd_i.dy2.*(IP*(rho_i'-rho_wl))) + ...
                                 + this.optsPhys.kBT*(rho_i(1) - rho_wl(1));                
                %dmuCheck(i) = -(kBT*(rho_i(1) - rho_wl(1)) - Int_1D*( (rho_i - rho_wl).*dVAdd.dy2))/(rhoLiq_eq-rhoGas_eq);
            end

            this.disjoiningPressureCheck = fB_AI;    
        end
    end
end