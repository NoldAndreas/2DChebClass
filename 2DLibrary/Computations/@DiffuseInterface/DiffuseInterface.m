classdef DiffuseInterface < handle
    
   properties (Access = public)              
       IC              
       optsNum,optsPhys                                    
   end
   
   
   methods (Access = public)          
        function this = DiffuseInterface(config)
            
            this.optsNum         = config.optsNum;
            this.optsPhys        = config.optsPhys;                                       
                        
        end
        
        function Preprocess(this)                        
            this.IC = InfCapillaryQuad(this.optsNum.PhysArea);    
            this.IC.ComputeAll(this.optsNum.PlotArea);               
            this.IC.SetUpBorders(this.optsNum.PhysArea.NBorder);            
        end                
        function rho = InitialGuessRho(this)
            PtsCart    = this.IC.GetCartPts();
            Cn         = this.optsPhys.Cn;
            
            rho        = - tanh(PtsCart.y1_kv/Cn);
        end        
        function theta = FindInterfaceAngle(this,rho)

            PtsCart      = this.IC.GetCartPts();
            fsolveOpts   = optimset('Display','off');                

            pt.y2_kv  =  min(PtsCart.y2_kv);
            y1CartStart = fsolve(@rhoX1,0,fsolveOpts);

            pt.y2_kv  =  max(PtsCart.y2_kv);
            y1CartEnd = fsolve(@rhoX1,0,fsolveOpts);                

            alpha = atan((y1CartStart-y1CartEnd)/...
                            (max(PtsCart.y2_kv)-  min(PtsCart.y2_kv)));

            theta = alpha + pi/2;

            function z = rhoX1(y1)
                pt.y1_kv = y1;
                IP       = this.IC.SubShapePtsCart(pt);
                z        = IP*rho;
            end    
        end        
                       
        %**************************************
        %Plot Functions
        %**************************************
        function PlotSeppecherSolution(this,D_B,theta,rho)
            
            UWall   = this.optsPhys.UWall;
            D_A     = this.optsPhys.D_A;
            PtsCart = this.IC.GetCartPts();
            y2Max   = this.optsNum.PhysArea.y2Max;            
                        
            u_flow = GetSeppecherSolutionCart(PtsCart,UWall,D_A,D_B,theta);                        
                        
            if(nargin >= 4)
                figure('Position',[0 0 1800 600],'Color','white');            
                subplot(1,2,1);  PlotU(this,u_flow);            
                subplot(1,2,2);  this.IC.doPlotFLine([-100,100],[y2Max,y2Max],rho.*u_flow(end/2+1:end),'CART'); 
            else
                figure('Position',[0 0 800 800],'Color','white');            
                PlotU(this,u_flow);            
            end
            title('Check accuracy of map');
        end        
        function PlotMu_and_U(this,mu,uv)            
            
            figure('Position',[0 0 1800 600],'color','white');
            subplot(1,2,1); this.IC.doPlots(mu);
            subplot(1,2,2); PlotU(this,uv);
        end        
        function PlotU(this,uv)            
            y2Max = this.optsNum.PhysArea.y2Max;
            
            y1Min = this.optsNum.PlotArea.y1Min;
            y1Max = this.optsNum.PlotArea.y1Max;
            
            y2L = (1:y2Max)';
            y1L = (y1Min:y1Max)';
            
            startPtsy1    = [y1Max*ones(size(y2L));...
                             y1Min*ones(size(y2L));...
                             y1L];
            startPtsy2    = [y2L;y2L;y2Max*ones(size(y1L))];
            this.IC.doPlotsStreamlines(uv,startPtsy1,startPtsy2); %IC.doPlotsFlux(u_flow)(mu);
        end
           
        %**********************************
        %External files:
        %**********************************
        rho = GetEquilibriumDensity(this,mu,theta,nParticles,rho)
        D_B = SetD_B(this,theta,rho,initialGuessDB)
        [mu,uv,A,b] = GetVelocityAndChemPot(this,rho,D_B,theta)
        [A,b] = ContMom_DiffuseInterfaceSingleFluid(this,rho)
        [A,b] = Div_FullStressTensor(this,rho)
        [A,b] = FullStressTensorIJ(this,rho,i,j)
   end
    
end