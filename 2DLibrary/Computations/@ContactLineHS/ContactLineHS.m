classdef ContactLineHS < DDFT_2D
    
    properties (Access = public)                                           
        %(1) 1D Results
        rho1D_lg,rho1D_wl,rho1D_wg,alpha_YCA  
        ST_1D  %1D surface tensions
        
        y1_SpectralLine        
        hII = [], hIII = [],hContour = []
        disjoiningPressure_II=[],disjoiningPressure_IV=[]
        y1_I = [], hI = []
        
        y2=[],Int_y2=[],DiffY2=[]
        
        IsolineInterfaceY2=[]        
               
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
             configuration.optsNum.PhysArea.shape = 'HalfSpace_FMT';
             this@DDFT_2D(configuration);
        end                
        
        %Preprocessing and testing of results
        function Preprocess(this)
            global dirDataOrg
            ChangeDirData([dirDataOrg filesep 'deg',num2str(this.optsNum.PhysArea.alpha_deg,3)]);
            Preprocess@DDFT_2D(this);     
            TestPreprocess(this);        
            
            %
            if(isfield(this.optsNum,'PlotAreaCart'))
                shapeSL = struct('yMin',this.optsNum.PlotAreaCart.y1Min,...
                                 'yMax',this.optsNum.PlotAreaCart.y1Max,...
                                 'N',100);                       
                this.y1_SpectralLine = SpectralLine(shapeSL);
                this.y1_SpectralLine.ComputeAll();
            end
            %
            ComputeST(this);
        end        
        function TestPreprocess(this)                          
                        
            R = this.IDC.R;             
             
            disp('*** Test preprocessed matrices ***');
                                    
            disp('*** Test convolution matrix ***');
            CheckMeanfieldConvolution(this);
            
            %(4) Test FMT Matrices            
            disp('*** Test FMT matrices ***');
            if(isfield(this.optsNum.FexNum,'Fex') && strcmp(this.optsNum.FexNum.Fex,'FMTRosenfeld_3DFluid'))
                CheckAverageDensities_Rosenfeld_3D(this.IDC,this.IntMatrFex);
            end
                        
            %(3) Test Integration
            disp('*** Test integration vector ***');
            a = 2;
            %(a) Integration in HalfSpace
            pts   = this.IDC.GetCartPts();
            r     = sqrt(pts.y1_kv.^2+(pts.y2_kv-R).^2);
            testF = exp(-(r/a).^2);    
            PrintErrorPos(this.IDC.Int*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);    

            %(b) Integration in Composed Half Space
            pts   = this.IDC.AD.GetCartPts();
            r     = sqrt(pts.y1_kv.^2+pts.y2_kv.^2);
            testF = exp(-(r/a).^2);
            Inth  = this.IDC.AD.ComputeIntegrationVector();
            PrintErrorPos(Inth*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);
            disp('*****');
                                    
        end         
                
        %1D surface tensions
        function [om,rho1D] = Compute1D(this,WLWGLG)
            rhoLiq_sat       = this.optsPhys.rhoLiq_sat;
            rhoGas_sat       = this.optsPhys.rhoGas_sat;
            optss            = this.optsPhys;   
            Fex_Num          = this.optsNum.FexNum;                        
            
            if(strcmp(WLWGLG,'WL'))
                optss.rho_iguess = this.optsPhys.rhoLiq_sat;
                [rho1D,params] = FMT_1D(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv);
                this.ST_1D.om_wallLiq = params.Fex;
                this.rho1D_wl         = rho1D;                
            elseif(strcmp(WLWGLG,'WG'))
                optss.rho_iguess = this.optsPhys.rhoGas_sat;
                [rho1D,params] = FMT_1D(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv);
                
                this.ST_1D.om_wallGas = params.Fex;
                this.rho1D_wg         = rho1D;                
            elseif(strcmp(WLWGLG,'LG'))
                Pts              = this.IDC.Pts;
                theta_CS         = this.IDC.alpha;  
            
                optss.rho_iguess = (rhoLiq_sat+rhoGas_sat)/2 + ...
                                  (rhoLiq_sat-rhoGas_sat)/2*tanh(Pts.y1*sin(theta_CS));
                [rho1D,params] = FMT_1D_Interface(this.IDC,this.IntMatrFex,optss,Fex_Num,this.IntMatrV2.Conv);                
                
                this.ST_1D.om_LiqGas  = params.Fex;
                this.rho1D_lg         = rho1D;
            end
            
            om  = params.Fex;                                  
            fprintf(['Surface tension ',WLWGLG,' = ',num2str(om),'\n']);
        end                
        function alpha_YCA = ComputeST(this)                        
            
            Compute1D(this,'WL');            
            Compute1D(this,'WG');
            Compute1D(this,'LG');                      

            alpha_YCA = ComputeContactAngle(this.ST_1D.om_wallGas,...
                            this.ST_1D.om_wallLiq,...
                            this.ST_1D.om_LiqGas);
            
            %****************************************
            this.alpha_YCA = alpha_YCA;
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
        function y0  = getInitialGuess(this,rho_ig)
                    
        	rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
            rhoGas_sat    = this.optsPhys.rhoGas_sat;
            kBT           = this.optsPhys.kBT;            
            
            p         = (this.rho1D_lg-rhoGas_sat)/(rhoLiq_sat-rhoGas_sat);    
            rho_ig    = kron(p,this.rho1D_wl) + kron(1-p,this.rho1D_wg);         
            y0        = kBT*log(rho_ig)+this.Vext;
        end
        
        %Plot functions        
        PlotContourResults(this,plain)
        function PlotDensityResult(this)                        
            rho           = GetRhoEq(this);
            
            figure('Color','white','Position',[0 0 1200 800]);
            this.IDC.doPlots(rho,'SC');
            zlabel('$\varrho$','Interpreter','Latex','fontsize',26);
            colormap(hsv);
            set(gca, 'CLim', [0, 1.0]);
            PlotArea = this.optsNum.PlotAreaCart;
            pbaspect([(PlotArea.y1Max-PlotArea.y1Min) (PlotArea.y2Max-PlotArea.y2Min) 5]);
            view([-10 5 3]);

            print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq],gcf);
            saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '.fig']);
            
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

            print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Interfaces'],gcf);
            saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Interfaces.fig']);

        end
    
        PlotDisjoiningPressures(this)
        PlotDensitySlices(this)
                
        %Adsorption Isotherm
        ComputeAdsorptionIsotherm(this,n,drying)
        FittingAdsorptionIsotherm(this,FT_Int,n)
        SumRule_AdsorptionIsotherm(this,ST_LG)
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
        function drying = IsDrying(this)
            drying = (min(this.AdsorptionIsotherm.FT) < 0);
        end
        
        %Compute functions 
        function ComputeEquilibrium(this)
            optsIn.maxComp_y2 = this.optsNum.maxComp_y2;
            miscIn.mark       = (this.IDC.GetCartPts.y2_kv <= ...
                                        this.optsNum.maxComp_y2);
            ComputeEquilibrium@DDFT_2D(this,[],optsIn,miscIn);
        end
        sol = Compute(this)        
        ComputeDynamics(this)        
        
        %Compute disjoining pressure        
        Compute_DisjoiningPressure_II(this,y1Int)
        Compute_DisjoiningPressure_IV(this)
        function SumRule_DisjoiningPressure_II(this,ST_LG)

            if(nargin < 2)
                ST_LG = this.ST_1D.om_LiqGas;
            end
            
            anaDP = this.disjoiningPressure_IV;

            %******************************************
            %Check sum rule Eq. (11) in Henderson (2005)    
            %Integrate along y1, and check sum rule
            % - Int(DisjoiningPressure(y),y=-inf..inf) = \gammaLV*sin(theta)            
            
            Int_y1  = this.y1_SpectralLine.Int;
            sinGamm = sin(this.alpha_YCA)*ST_LG;
            err       = (Int_y1*anaDP + sinGamm);
            errPC     = err/sinGamm*100;            
            disp(['Integral expression: ',num2str(Int_y1*anaDP),' +/- ',num2str(err) , ' or relative: ',num2str(errPC),' percent']);

            h              = (-Int_y1*anaDP)/ST_LG;
            estTheta       = asin(h)*180/pi;
            if(IsDrying(this))
                estTheta = 180 - estTheta;
            end
            error_estTheta = 180/pi*sum(Int_y1)*max(abs(anaDP(1)),abs(anaDP(end)))/ST_LG/cos(estTheta*pi/180);%*(1/sqrt(1+h^2));
            disp(['Theta from Sum rule = ',num2str(estTheta),' [deg] +/- ',num2str(err/ST_LG/cos(estTheta*pi/180)*180/pi),' [deg]']);                            
            disp(['Error from numerical estimate : ',num2str(error_estTheta),' [deg]']);
            %disp(['Error from difference to integral estimate : ',num2str(err/ST_LG/cos(estTheta*pi/180)*180/pi),' [deg]']);

            %PrintErrorPos(180/pi*(estTheta-this.alpha_YCA),'Estimated contact angle through sum rule integrating disjoining pressure [percent]');
        end        
        
        %Compute height profiles
        Compute_hI(this)
        Compute_hII(this)
        Compute_hIII(this)
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
            
            hS           = max(this.hI);
            h0           = min(this.hIII);

            fsolveOpts   = optimset('Display','off');            
            f            = this.hIII-h0;
            [~,j]        = max(this.hIII);
            [DeltaY1_III,~,exitflag]  = fsolve(@fX,this.y1_SpectralLine.Pts.y(j),fsolveOpts);            
            if(exitflag < 1)
                cprintf('*r','ComputeDeltaFit: Fitting hIII vs hI: no solution found');
            end
            
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
        
        %to delete                   
        [f,y1]  = PostProcess(this,y1Int)                                
                                
        I = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)        
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);
        
        PlotDisjoiningPressureAnalysis(this)    
        PlotInterfaceAnalysisY1(this)
        [y2,theta] = PlotInterfaceAnalysisY2(this,yInt)                
        PlotDensitySlicesMovie(this);       

    end
end