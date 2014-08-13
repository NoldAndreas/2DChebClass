classdef ContactLineHS < DDFT_2D
    
    properties (Access = public)                   
        configName
                
        %(1) 1D Results
        rho1D_lg,rho1D_wl,rho1D_wg,alpha_YCA  
        ST_1D  %1D surface tensions
        
        %(2) 2D Results
        FilenameEq
        
        y2=[],Int_y2=[],DiffY2=[]
        IsolineInterfaceY2=[]
        y1=[],Int_y1=[],DiffY1=[],DiffYY1=[]
        filmThickness = [],IsolineInterface = []
        ell_2DDisjoiningPressure=[]
               
        %structure with information of the adsorption isotherm        
        AdsorptionIsotherm 
        
        grandPot,grandPot2
        disjoiningPressure,disjoiningPressureCheck
        
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
            ComputeST(this);
        end        
        function TestPreprocess(this)                          
            
             M = this.IDC.M;
             R = this.IDC.R;             
                                    
             %(3.1) Test Convolution at infinity
             fMF       = str2func(this.optsPhys.V2.V2DV2);
             [h1,h2,a] = fMF(1,this.optsPhys.V2);
             marky2Inf = (this.IDC.Pts.y2_kv == inf);
             PrintErrorPos(this.IntMatrV2.Conv(marky2Inf,:)*ones(M,1)- 2*a,'Error for convolution at y2 = infinity',this.IDC.Pts.y1_kv(marky2Inf));          

             %Convolution profile
             if(strcmp(this.optsPhys.V2.V2DV2,'Phi2DLongRange'))
                 y0R = PtsCart.y2_kv-R;
                 h = this.optsPhys.V2.epsilon*(- pi^2/2 + ...
                                          + pi*atan(-y0R)+...
                                          - pi*(y0R)./(1+y0R.^2));
                 PrintErrorPos(h-this.IntMatrV2.Conv*ones(M,1),'Error of Phi2DLongRange*1',this.IDC.GetCartPts);
             elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHenderson_2D'))                 
                 conv = this.IntMatrV2.Conv(this.IDC.Pts.y1_kv==inf,:);
                 y2_h = this.IDC.GetCartPts.y2_kv(this.IDC.Pts.y1_kv==inf) - R;                          
                 Psi(y2_h < 1)  = -16/9*pi +6/5*pi*y2_h(y2_h < 1);         
                 Psi(y2_h >= 1) = 4*pi*(1./(45*y2_h(y2_h >= 1).^9) - 1./(6*y2_h(y2_h >= 1).^3));
                 check          = 2*a-this.optsPhys.V2.epsilon*Psi;
                 
                 PrintErrorPos(conv*ones(M,1) - check','Error for convolution at y1 = infinity',y2_h);                 
             end     
            
            %(4) Test FMT Matrices            
            if(isfield(this.optsNum.FexNum,'Fex') && strcmp(this.optsNum.FexNum.Fex,'FMTRosenfeld_3DFluid'))
                CheckAverageDensities_Rosenfeld_3D(this.IDC,this.IntMatrFex);
            end
            
            %(3) Test Integration
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
        function InitInterpolation(this,Nodefault,PlotArea)
            
            if((nargin == 1) || ~Nodefault)
                if(abs(this.optsNum.PhysArea.alpha_deg - 90) < 20)
                    PlotArea = struct('y1Min',-5,'y1Max',5,'y2Min',0.5,'y2Max',10,'N1',80,'N2',80);
                    %PlotArea = struct('y1Min',-8,'y1Max',8,...
                    %                  'y2Min',0.5,'y2Max',this.optsNum.maxComp_y2+5,...
                    %                  'N1',80,'N2',80);
                elseif(this.optsNum.PhysArea.alpha_deg <= 70)
                    if(isempty(this.optsNum.maxComp_y2))
                        y2MaxPlot = 10;
                    else
                        y2MaxPlot = this.optsNum.maxComp_y2+5;
                    end
                    PlotArea = struct('y1Min',-2,...
                                      'y1Max',y2MaxPlot/tan(this.optsNum.PhysArea.alpha_deg*pi/180),...
                                      'y2Min',0.5,'y2Max',y2MaxPlot,...
                                      'N1',100,'N2',100);        
                end    
                this.optsNum.PlotArea = PlotArea;
                this.IDC.InterpolationPlotCart(PlotArea,true);
            else
                if(nargin > 2)
                    if(~isfield(this.optsNum,'PlotArea') || ~isequal(this.optsNum.PlotArea,PlotArea))
                        this.optsNum.PlotArea = PlotArea;
                        this.HS.InterpolationPlotCart(this.optsNum.PlotArea,true);
                    end
                else
                    this.HS.InterpolationPlotCart(this.optsNum.PlotArea,true);                    
                end                
            end
                        
        end                
        function SetV1(this,epw)
            this.optsPhys.V1.epsilon_w = epw;
            PtsCart                    = this.HS.GetCartPts();
            this.VAdd                  = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,this.optsPhys.V1);
        end       
        function InitAnalysisGrid(this,y1Int,y2Int)
            if(~isempty(y1Int))
                [this.y1,this.Int_y1,this.DiffY1,this.DiffYY1] = InitAnalysisGridY(this,y1Int,100,'CHEB');        
            end
            if(~isempty(y2Int))
                [this.y2,this.Int_y2,this.DiffY2] = InitAnalysisGridY(this,y2Int,100);                      
            end
        end        
        
        %Plot functions
        [fContour] =  PlotEquilibriumResults(this,bounds1,bounds2,plain,saveFigs)
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
        sol = Compute(this)        
        ComputeDynamics(this)        
                        
        %Postprocess functions        
        [f,y1]   = PostProcess(this,y1Int)
        PostProcess_2DDisjoiningPressure(this,y1Int)
        PostProcess_FilmThickness(this,y1Int)
        [f1,f2]           = Post_HFrom2DDisjoiningPressure(this,f1)
        [y2Cart]          = ComputeInterfaceContour(this,level)
        [y1Cart]          = ComputeInterfaceContourY2(this,level,y2)
        [CA_measured,err] = MeasureContactAngle(this,type,yInt)
        
                        
        I = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)
        SumRule_DisjoiningPotential(this,ST_LG)
        
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);
    end
end