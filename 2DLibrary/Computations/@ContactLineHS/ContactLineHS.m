classdef ContactLineHS < DDFT_2D
    
    properties (Access = public)                   
        configName
                
        %(1) 1D Results
        rho1D_lg,rho1D_wl,rho1D_wg,alpha_YCA  
        ST_1D  %1D surface tensions
        
        y1_SpectralLine        
        hII = [], hIII = [],hContour = []
        disjoiningPressure_II
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
            
             global dirData
             if(nargin > 0)                
                this.configName    = SaveConfig(configuration,'Configurations');        
             else
                [configIn,DataFolder] = uigetfile([dirData filesep 'Configurations' filesep '*.mat'],['Select Config File']);
                load([DataFolder,configIn]);
                disp(['Loading configuration ',[DataFolder,configIn],' ...']);
                this.configName = configIn(1:end-4);        
             end   
                                                    
        end                
        
        %Preprocessing and testing of results
        function Preprocess(this)
            global dirDataOrg
            ChangeDirData([dirDataOrg filesep 'deg',num2str(this.optsNum.PhysArea.alpha_deg,3)]);
            Preprocess@DDFT_2D(this);     
            TestPreprocess(this);        
            
            %
            shapeSL = struct('yMin',this.optsNum.PlotAreaCart.y1Min,...
                             'yMax',this.optsNum.PlotAreaCart.y1Max,...
                             'N',100);                       
            this.y1_SpectralLine = SpectralLine(shapeSL);
            this.y1_SpectralLine.ComputeAll();
            %
            ComputeST(this);
        end        
        function TestPreprocess(this)                          
            
             M = this.IDC.M;
             R = this.IDC.R;             
             
             disp('*** Test preprocessed matrices ***');
                                    
             disp('*** Test convolution matrix ***');
             %Test Convolution at infinity             
             fMF       = str2func(this.optsPhys.V2.V2DV2);
             [h1,h2,a] = fMF(1,this.optsPhys.V2);
             marky2Inf = (this.IDC.Pts.y2_kv == inf);
             PrintErrorPos(this.IntMatrV2.Conv(marky2Inf,:)*ones(M,1)- 2*a,'convolution at y2 = infinity',this.IDC.Pts.y1_kv(marky2Inf));          

             %Convolution profile
             if(strcmp(this.optsPhys.V2.V2DV2,'Phi2DLongRange'))
                 y0R = PtsCart.y2_kv-R;
                 h = this.optsPhys.V2.epsilon*(- pi^2/2 + ...
                                          + pi*atan(-y0R)+...
                                          - pi*(y0R)./(1+y0R.^2));
                 PrintErrorPos(h-this.IntMatrV2.Conv*ones(M,1),'Phi2DLongRange*1',this.IDC.GetCartPts);
             elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHenderson_2D'))                 
                 conv = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
                 y2_h = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - R;                          
                 Psi(y2_h < 1)  = -16/9*pi +6/5*pi*y2_h(y2_h < 1);         
                 Psi(y2_h >= 1) = 4*pi*(1./(45*y2_h(y2_h >= 1).^9) - 1./(6*y2_h(y2_h >= 1).^3));
                 check          = 2*a-this.optsPhys.V2.epsilon*Psi;
                 
                 PrintErrorPos(conv*ones(M,1) - check','convolution at y1 = infinity',y2_h);                 
             end     
            
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
        function InitAnalysisGrid(this,y2Int)
            

            if(~isempty(y2Int))
                [this.y2,this.Int_y2,this.DiffY2] = InitAnalysisGridY(this,y2Int,100);                      
            end
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
        [fContour] =  PlotEquilibriumResults(this,plain,saveFigs)
        
        PlotDisjoiningPressures(this)
        
        PlotDisjoiningPressureAnalysis(this)    
        PlotInterfaceAnalysisY1(this)
        [y2,theta] = PlotInterfaceAnalysisY2(this,yInt)        
        PlotDensitySlices(this);
        PlotDensitySlicesMovie(this);       
        
        %Adsorption Isotherm
        ComputeAdsorptionIsotherm(this,n,drying)
        FittingAdsorptionIsotherm(this,FT_Int,n)
        SumRule_AdsorptionIsotherm(this,ST_LG)
        
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
        
        %Compute height profiles
        Compute_hI(this)
        Compute_hII(this)
        Compute_hIII(this)
        [y2Cart] = Compute_hContour(this,level)        
        
        %Postprocess
        [CA_measured,err] = MeasureContactAngle(this,type,yInt)
        [y1Cart]          = ComputeInterfaceContourY2(this,level,y2)        
        
        %to delete      
        Post_HFrom2DDisjoiningPressure(this,f1)        
        [f,y1]  = PostProcess(this,y1Int)                                
        
                        
        I = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)
        SumRule_DisjoiningPotential(this,ST_LG)        
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);
    end
end