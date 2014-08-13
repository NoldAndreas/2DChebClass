classdef ContactLine < handle
    
    properties (Access = public)                   
        optsNum
        optsPhys
        configName
        
        HS
        Conv
        IntMatrFex
        Vext,Vext_grad,VAdd
        
        %Computational Results:
        %(1) 1D Results
        rho1D_lg,rho1D_wl,rho1D_wg,alpha_YCA  
        ST_1D  %1D surface tensions
        
        %(2) 2D Results
        rho_eq,FilenameEq
        
        
        y2=[],Int_y2=[],DiffY2=[]
        IsolineInterfaceY2=[]
        y1=[],Int_y1=[],DiffY1=[],DiffYY1=[]
        filmThickness = [],IsolineInterface = []
        ell_2DDisjoiningPressure=[]
       
        
        AdsorptionIsotherm_FT,AdsorptionIsotherm_Mu,AdsorptionIsotherm_rho,AdsorptionIsotherm_OmEx,AdsorptionIsotherm_dmuCheck,AdsorptionIsotherm_Coeff
        AdsorptionIsotherm_Pts
        grandPot,grandPot2
        disjoiningPressure,disjoiningPressureCheck
        
        wIntLoc,wIntHS,f_loc,f_hs
        CA_deg_measured,CA_deg_measured_error
        
        %***************
        %Other Parameters
        saveFigs = true;
    end
    
    methods (Access = public)          
        function this = ContactLine(configuration)
             global dirData
             if(nargin > 0)                
                this.configName    = SaveConfig(configuration,'Configurations');        
             else
                [configIn,DataFolder] = uigetfile([dirData filesep 'Configurations' filesep '*.mat'],['Select Config File']);
                load([DataFolder,configIn]);
                disp(['Loading configuration ',[DataFolder,configIn],' ...']);
                this.configName = configIn(1:end-4);        
             end   
             
             this.optsNum         = configuration.optsNum;
             this.optsPhys        = configuration.optsPhys;                          
             this.optsPhys.HSBulk = (['FexBulk_',this.optsNum.FexNum.Fex]);      
             
             
             theta_CS = this.optsNum.PhysArea.alpha_deg*pi/180;
             %y1Min   = 0;
             %y1Max   = this.optsNum.maxComp_y2/tan(theta_CS);%PlotArea.y1Max+3;      
     
             %m       = 100;
             %this.y1 = y1Min + (0:(m-1))'/(m-1)*(y1Max-y1Min);             
        end                
        function Preprocess(this)             
            %************************************************
            %****************  Preprocess  ****************
            %************************************************                               
            global MinimalOutput
            GotoSubDir(this);
                        
            N1 = this.optsNum.PhysArea.N(1);
            N2 = this.optsNum.PhysArea.N(2);            
            R   = this.optsPhys.sigmaS/2;     

            %(1) Thermodynamic Values
            [rhoGas_sat,rhoLiq_sat,mu_sat] = ...
                BulkSatValues(this.optsPhys,[0.01;0.6;-2],~MinimalOutput);
            GetCriticalPoint(this.optsPhys);
           % BulkPhaseDiagram(this.optsPhys);

            this.optsPhys.mu_sat      = mu_sat;
            this.optsPhys.rhoGas_sat  = rhoGas_sat;
            this.optsPhys.rhoLiq_sat  = rhoLiq_sat;

            %(2) Numerical Integration, Differentiation
            optsHS             = this.optsNum.PhysArea;
            optsHS.alpha       = this.optsNum.PhysArea.alpha_deg*pi/180;
            this.HS            = HalfSpace_FMT(optsHS,diag(this.optsPhys.sigmaS)/2);
            [Pts,Diff,Int,Ind] = this.HS.ComputeAll();
            PtsCart            = this.HS.GetCartPts();

            %(3) Numerical Convolution
            opts.V2                  = this.optsPhys.V2;
            opts.nSpecies            = this.optsPhys.nSpecies;
            opts.optsNum.PhysArea.N  = this.optsNum.PhysArea.N;
            opts.optsNum.PhysArea.L1 = this.optsNum.PhysArea.L1;
            opts.optsNum.PhysArea.L2 = this.optsNum.PhysArea.L2;
            opts.optsNum.PhysArea.alpha_deg = this.optsNum.PhysArea.alpha_deg;
            opts.optsNum.PhysArea.Conv      = this.optsNum.PhysArea.Conv;
            opts.Comments           = ['ConfigFile: ',this.configName,'.txt'];
            convStruct              = DataStorage(['HalfSpace_FMT' filesep 'FexMatrices_SplitDisk'],@FexMatrices_SplitDisk,opts,this.HS);
            this.Conv               = convStruct.Conv;

             %(3.1) Test Convolution    
             %at infinity
             fMF       = str2func(this.optsPhys.V2.V2DV2);
             [h1,h2,a] = fMF(1,this.optsPhys.V2);
             marky2Inf = (this.HS.Pts.y2_kv == inf);
             PrintErrorPos(this.Conv(marky2Inf,:)*ones(N1*N2,1)- 2*a,'Error for convolution at y2 = infinity',this.HS.Pts.y1_kv(marky2Inf));          

             %Convolution profile
             if(strcmp(this.optsPhys.V2.V2DV2,'Phi2DLongRange'))
                 y0R = PtsCart.y2_kv-R;
                 h = this.optsPhys.V2.epsilon*(- pi^2/2 + ...
                                          + pi*atan(-y0R)+...
                                          - pi*(y0R)./(1+y0R.^2));
                 PrintErrorPos(h-this.Conv*ones(N1*N2,1),'Error of Phi2DLongRange*1',PtsCart);
             elseif(strcmp(this.optsPhys.V2.V2DV2,'BarkerHenderson_2D'))                 
                 conv = this.Conv(Pts.y1_kv==inf,:);
                 y2_h = PtsCart.y2_kv(Pts.y1_kv==inf) - R;                          
                 Psi(y2_h < 1)  = -16/9*pi +6/5*pi*y2_h(y2_h < 1);         
                 Psi(y2_h >= 1) = 4*pi*(1./(45*y2_h(y2_h >= 1).^9) - 1./(6*y2_h(y2_h >= 1).^3));
                 check          = 2*a-this.optsPhys.V2.epsilon*Psi;
                 
                 PrintErrorPos(conv*ones(N1*N2,1) - check','Error for convolution at y1 = infinity',y2_h);                 
             end     
            %*******************************************         
            %(4) FMT Matrices            
            ComputeHardSphereMatrices(this);

             %(5) External Potential           
            [this.Vext,this.Vext_grad]  = getVBackDVBack(PtsCart.y1_kv,PtsCart.y2_kv,this.optsPhys.V1);                  
            this.VAdd                   = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,this.optsPhys.V1);

            %(3) Test Integration
            a = 2;
            %(a) Integration in HalfSpace
            pts   = this.HS.GetCartPts();
            r     = sqrt(pts.y1_kv.^2+(pts.y2_kv-R).^2);
            testF = exp(-(r/a).^2);    
            PrintErrorPos(this.HS.Int*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);    

            %(b) Integration in Composed Half Space
            pts   = this.HS.AD.GetCartPts();
            r     = sqrt(pts.y1_kv.^2+pts.y2_kv.^2);
            testF = exp(-(r/a).^2);
            Inth  = this.HS.AD.ComputeIntegrationVector();
            PrintErrorPos(Inth*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);
            
            %*******************************************         
            ComputeST(this);
        end         
        
        ComputeHardSphereMatrices(this)        
        ComputeAdsorptionIsotherm(this,n,drying)    
        
        function [om,rho1D] = Compute1D(this,plotBool,WLWGLG)
            rhoLiq_sat       = this.optsPhys.rhoLiq_sat;
            rhoGas_sat       = this.optsPhys.rhoGas_sat;
            optss            = this.optsPhys;   
            Fex_Num          = this.optsNum.FexNum;                        
            
            if(strcmp(WLWGLG,'WL'))
                optss.rho_iguess = this.optsPhys.rhoLiq_sat;
                [rho1D,parms] = FMT_1D(this.HS,this.IntMatrFex,optss,Fex_Num,this.Conv,plotBool);
            elseif(strcmp(WLWGLG,'WG'))
                optss.rho_iguess = this.optsPhys.rhoGas_sat;
                [rho1D,parms] = FMT_1D(this.HS,this.IntMatrFex,optss,Fex_Num,this.Conv,plotBool);
            elseif(strcmp(WLWGLG,'LG'))
                Pts              = this.HS.Pts;
                theta_CS         = this.optsNum.PhysArea.alpha_deg*pi/180; 	            
            
                optss.rho_iguess = (rhoLiq_sat+rhoGas_sat)/2 + ...
                                  (rhoLiq_sat-rhoGas_sat)/2*tanh((Pts.y1-this.optsNum.y1Shift)*sin(theta_CS));
                [rho1D,parms] = FMT_1D_Interface(this.HS,this.IntMatrFex,optss,Fex_Num,this.Conv,plotBool,this.optsNum.y1Shift);                
            end
            om  = parms.Fex;           
        end        
        function [alpha_YCA,ST_1D,rho1D_lg,rho1D_wl,rho1D_wg] = ComputeST(this,plotBool)
            
            global MinimalOutput                        
            if(nargin == 1)
                plotBool = ~MinimalOutput;
            end            
            GotoSubDir(this);
            %**********
            
            [ST_1D.om_wallLiq,rho1D_wl] = Compute1D(this,plotBool,'WL');            
            [ST_1D.om_wallGas,rho1D_wg] = Compute1D(this,plotBool,'WG');
            [ST_1D.om_LiqGas,rho1D_lg]  = Compute1D(this,plotBool,'LG');                      
 
            if(~MinimalOutput)
                fprintf(['Surface tension (Liq/Gas) = ',num2str(ST_1D.om_LiqGas),'\n']);
                fprintf(['Surface tension (Wall/Liq) = ',num2str(ST_1D.om_wallLiq),'\n']);
                fprintf(['Surface tension (Wall/Gas) = ',num2str(ST_1D.om_wallGas),'\n']);
            end

            alpha_YCA = ComputeContactAngle(ST_1D.om_wallGas,ST_1D.om_wallLiq,ST_1D.om_LiqGas);
            
            %****************************************
            this.ST_1D = ST_1D;
            this.rho1D_wg = rho1D_wg;
            this.rho1D_wl = rho1D_wl;
            this.rho1D_lg = rho1D_lg;
            this.alpha_YCA = alpha_YCA;
        end        
        function GotoSubDir(this)
            global dirDataOrg
            ChangeDirData([dirDataOrg filesep 'deg',num2str(this.optsNum.PhysArea.alpha_deg,3)]);
        end        
        function InitInterpolation(this,PlotArea)
            
            if((nargin == 1) || ( (ischar(PlotArea) == 2) && strcmp(PlotArea,'newDefault')))
                if(abs(this.optsNum.PhysArea.alpha_deg - 90) < 20)
                    PlotArea = struct('y1Min',-5,'y1Max',5,'y2Min',0.5,'y2Max',10,'N1',80,'N2',80);
                    %PlotArea = struct('y1Min',-8,'y1Max',8,...
                    %                  'y2Min',0.5,'y2Max',this.optsNum.maxComp_y2+5,...
                    %                  'N1',80,'N2',80);
                elseif((this.optsNum.PhysArea.alpha_deg -90) > 20)
                    if(isempty(this.optsNum.maxComp_y2))
                        y2MaxPlot = 10;
                    else
                        y2MaxPlot = this.optsNum.maxComp_y2+5;
                    end                    
                    PlotArea = struct('y2Min',0.5,'y2Max',y2MaxPlot,...
                                      'N1',100,'N2',100);        
                    if(this.optsNum.PhysArea.alpha_deg < 90)
                        PlotArea.y1Min = -2;
                        PlotArea.y1Max = y2MaxPlot/tan(this.optsNum.PhysArea.alpha_deg*pi/180);
                    else
                        PlotArea.y1Min = y2MaxPlot/tan(this.optsNum.PhysArea.alpha_deg*pi/180);
                        PlotArea.y1Max = 10;
                    end
                                  
                end    
                this.optsNum.PlotArea = PlotArea;
                this.HS.InterpolationPlotCart(PlotArea,true);
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
        
        sol = Compute(this)
        ComputeEquilibrium(this,redo)
        ComputeDynamics(this)        
                        
        %Postprocess functions
        InitAnalysisGrid(this,y1Int,y2Int)
        [f,y1]   = PostProcess(this,y1Int)
        PostProcess_2DDisjoiningPressure(this,y1Int)
        PostProcess_FilmThickness(this,y1Int)
        [f1,f2]           = Post_HFrom2DDisjoiningPressure(this,f1)
        [y2Cart]          = ComputeInterfaceContour(this,level)
        [y1Cart]          = ComputeInterfaceContourY2(this,level,y2)
        [CA_measured,err] = MeasureContactAngle(this,type,yInt)
        FittingAdsorptionIsotherm(this,FT_Int,n)
        
        [fContour] =  PlotEquilibriumResults(this,bounds1,bounds2,plain,saveFigs)
        PlotDisjoiningPressureAnalysis(this)    
        PlotInterfaceAnalysisY1(this)
        [y2,theta] = PlotInterfaceAnalysisY2(this,yInt)
        
        PlotDensitySlices(this);
        PlotDensitySlicesMovie(this);
                
        I = doIntNormalLine(this,y2Max,y1,f_loc,f_hs)
        SumRule_DisjoiningPotential(this,ST_LG)
        
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);
    end
end