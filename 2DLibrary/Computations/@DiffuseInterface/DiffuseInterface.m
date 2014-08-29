classdef DiffuseInterface < handle
    
   properties (Access = public)              
       IC              
       optsNum,optsPhys    
       IntSubArea
       PlotBCShape
       
       rho = [],uv  = [],theta
       %surfaceTension = 4/3;
   end
   
   
   methods (Access = public)          
        function this = DiffuseInterface(config)             
            this.optsNum         = config.optsNum;
            this.optsPhys        = config.optsPhys;                                                               
        end       
        function Preprocess(this)                        
            this.IC = InfCapillaryQuad(this.optsNum.PhysArea);    
            
            this.IC.ComputeIndices();
            this.IC.ComputeDifferentiationMatrix();
            this.IC.ComputeIntegrationVector();
            this.IC.InterpolationPlot(this.optsNum.PlotArea,true);                          
            
            Sel = {'Dy1' ;'DDy1' ; 'Dy2'; 'DDy2';...
                   'DDDy2';'DDDy1';...
                   'Dy1Dy2'; 'DDy1Dy2'; 'Dy1DDy2';...
                   'Lap' ;'grad' ;'div'; ...
                   'gradLap' ;'gradDiv'; 'LapVec';'LapDiv'};
            this.IC.ComputeDifferentiationMatrix(Sel);
            
            this.IC.SetUpBorders(this.optsNum.PhysArea.NBorder);            

            bxArea          = this.optsNum.PlotArea;
            bxArea.N        = [100,100];
            BX              = Box(bxArea);
            IntBx           = BX.ComputeIntegrationVector();
            this.IntSubArea = IntBx*this.IC.SubShapePts(BX.GetCartPts());
            
            shapeBC         = this.optsNum.PhysArea;
            shapeBC.y2Min   = shapeBC.y2Max;
            shapeBC.y2Max   = shapeBC.y2Max + 1;
            this.PlotBCShape  = InfCapillaryQuad(shapeBC);
            
            plotBC          = this.optsNum.PlotArea;
            plotBC.y2Min    = plotBC.y2Max;
            plotBC.y2Max    = plotBC.y2Max+1;
            plotBC.N1 = 10; plotBC.N2 = 2;
            this.PlotBCShape.ComputeAll(plotBC);
        end        
        function rho = InitialGuessRho(this)
            PtsCart    = this.IC.GetCartPts();
            Cn         = this.optsPhys.Cn;
            
            rho        = tanh(PtsCart.y1_kv/Cn);
        end                      
                    
        SolveMovingContactLine(this,maxIterations)
        [rho,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,rho,findTheta)
        [rho,muDelta] = GetEquilibriumDensityR(this,mu,theta,nParticles,rho,ptC)                  
        
        [mu,uv,A,b,a] = GetVelocityAndChemPot(this,rho,theta)
        [A,b]       = ContMom_DiffuseInterfaceSingleFluid(this,rho)
        [rho,uv]    = SolveFull(this,ic)
        
        deltaX = GetDeltaX(this,rho,theta)              
        [uvBound_Corr,a] = CorrectVelocityProfile(this,theta,rho)
                         
        function [bulkError,bulkAverageError] = DisplayFullError(this,rho,uv)            
            M        = this.IC.M;
            
            [Af,bf]  = ContMom_DiffuseInterfaceSingleFluid(this,rho);                         
            mu_s     = GetMu(this,rho);
            
            error    = Af*[mu_s;uv] - bf;
            PrintErrorPos(error(1:M),'continuity equation',this.IC.Pts);            
            PrintErrorPos(error(1+M:2*M),'y1-momentum equation',this.IC.Pts);            
            PrintErrorPos(error(2*M+1:end),'y2-momentum equation',this.IC.Pts);                        
            
            PrintErrorPos(error(repmat(~this.IC.Ind.bound,2,1)),'bulk continuity and momentum equations');               
            
            bulkError        = max(abs(error(repmat(~this.IC.Ind.bound,2,1))));
            bulkAverageError = mean(abs(error(repmat(~this.IC.Ind.bound,2,1))));
        end        
        function spt = FindStagnationPoint(this)
            uv = this.uv;
            fsolveOpts = optimset('Display','off');
            [spt,~,flag] = fsolve(@ValueAtXY,[-2,2],fsolveOpts);                        
            if(flag < 1)
                disp('No stagnation point found');
                spt = [];
            else
                disp(['Stagnation point found at (',num2str(spt(1)),',',num2str(spt(2)),')']);
            end
            function z = ValueAtXY(xy)
                pt.y1_kv = xy(1);
                pt.y2_kv = xy(2);
                IP  = this.IC.SubShapePtsCart(pt);
                z   = [IP*uv(1:end/2);IP*uv(1+end/2:end)];
            end
        end
        
        %Auxiliary Cahn-Hilliard
        [A,b]       = Div_FullStressTensor(this,rho)
        [A,b]       = FullStressTensorIJ(this,rho,i,j)   
        function mu_s = GetMu(this,rho)
            if(nargin == 1)
                rho = this.rho;
            end
            Cn       = this.optsPhys.Cn;
            mu_s     = DoublewellPotential(rho,Cn) - Cn*(this.IC.Diff.Lap*rho);               
        end        
        function p = GetPressure_from_ChemPotential(this,mu,rho_ig)
            Cn     = this.optsPhys.Cn;
            rho_m  = this.optsPhys.rho_m;
            %for bulk
            % 0 = W' - m
            % p = - W + mu*(rho + rho_m)
            
            fsolveOpts = optimset('Display','off');
            [rho,~,exitflag]   = fsolve(@f,rho_ig,fsolveOpts);
            if(exitflag < 1)
                cprintf('*r','No solution for density for given chemical potential found');
            end
            
            [~,W] = DoublewellPotential(rho,Cn);
            p     = - W + mu*(rho+rho_m);
            
            function y = f(rho)
                y = DoublewellPotential(rho,Cn) - mu;
            end
            
        end                   
        
        %Plotting                           
        function PlotResultsMu(this,mu,uv) 
            if(nargin == 1)
                uv  = this.uv;
                mu  = GetMu(this,this.rho);            
            end
            
            figure('Position',[0 0 800 600],'color','white');
            this.IC.doPlots(mu,'contour');  
            PlotU(this,uv); hold on; 
        end
        function PlotResultsRho(this,uv,rho,theta)
            if(nargin == 1)
                uv  = this.uv;
                rho = this.rho;
                theta = this.theta;
            end
            figure('Position',[0 0 800 600],'color','white');
            PlotU(this,uv); hold on; 
            this.IC.doPlots(rho,'contour');     
            sp = FindStagnationPoint(this);hold on;
            plot(sp(1),sp(2),'or','MarkerFaceColor','r','MarkerSize',10); 
            
            this.PlotSeppecherSolution(theta,rho);
        end        
        function PlotSeppecherSolution(this,theta,rho)
            if(nargin == 1)                
                rho = this.rho;            
                theta = this.theta;
            end
            UWall   = this.optsPhys.UWall;
%            D_A     = this.optsPhys.D_A;                        
%             u_flow = GetSeppecherSolutionCart(this.PlotBCShape.GetCartPts,...
%                                               UWall,D_A,D_B,theta);                                                                        
%             this.PlotBCShape.doPlots(u_flow,'flux',struct('reshape',false,'linecolor','m'));
            
            PtsCart = this.IC.GetCartPts(); 
            if(nargin > 2)
                deltaX = GetDeltaX(this,rho,theta);             
                PtsCart.y1_kv = PtsCart.y1_kv - deltaX;
            else
                deltaX = 0;
            end             
             
             %y2Max   = this.optsNum.PhysArea.y2Max;            
                         
             u_flow = GetSeppecherSolutionCart(PtsCart,UWall,0,0,theta);     
%                         
%             if(nargin >= 4)
%                 figure('Position',[0 0 1800 600],'Color','white');            
%                 subplot(1,2,1);  
            PlotU(this,u_flow);            hold on;                   
            
            y2Max = this.optsNum.PhysArea.y2Max;
            plot([deltaX (deltaX+y2Max/tan(theta))],[0 y2Max],'k--','linewidth',2.5);
%                 subplot(1,2,2);  this.IC.doPlotFLine([-100,100],[y2Max,y2Max],rho.*u_flow(end/2+1:end),'CART'); 
%             else
%                 figure('Position',[0 0 800 800],'Color','white');            
%                 PlotU(this,u_flow);            
%             end
%             title('Check accuracy of map');
        end          
        function SavePlotResults(this,uv,rho,theta,mu)
            if(nargin == 1)
                uv  = this.uv;
                rho = this.rho;
                mu  = GetMu(this,rho);
                theta = this.theta;
            end
            
            global dirData
            
            filename = getTimeStr();
            
            PlotResultsRho(this,uv,rho,theta);
            print2eps([dirData filesep 'Density_' filename ],gcf);
            saveas(gcf,[dirData filesep 'Density_' filename '.fig']);
            
            PlotResultsMu(this,mu,uv);
            print2eps([dirData filesep 'ChemPot' filename],gcf);
            saveas(gcf,[dirData filesep 'ChemPot' filename '.fig']);
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
        
        %Old
        D_B         = SetD_B(this,theta,rho,initialGuessDB)
        function Interp = ResetOrigin(this,rho)
            
            pt.y2_kv  =  0;
            DeltaY1   = fsolve(@rhoX1,-1);%,fsolveOpts);
            
            ptsCartShift       = this.IC.GetCartPts();
            ptsCartShift.y1_kv = ptsCartShift.y1_kv + DeltaY1;
            
            Interp = this.IC.SubShapePtsCart(ptsCartShift);                       
            
            function z = rhoX1(y1)
                pt.y1_kv = y1;
                IP       = this.IC.SubShapePtsCart(pt);
                z        = IP*rho;
            end    
        end
        function theta = FindInterfaceAngle(this,rho)

            y2M = 0; y2P = 10;
          %  fsolveOpts   = optimset('Display','off');                

            pt.y2_kv  =  y2M;
            y1CartStart = fsolve(@rhoX1,-1);%,fsolveOpts);

            pt.y2_kv  = y2P;
            y1CartEnd = fsolve(@rhoX1,2);%,fsolveOpts);                

            alpha = atan((y1CartStart-y1CartEnd)/(y2P-  y2M));
            theta = alpha + pi/2;
            
            disp(['Contact angle: ',num2str(theta*180/pi),'[deg].']);

            function z = rhoX1(y1)
                pt.y1_kv = y1;
                IP       = this.IC.SubShapePtsCart(pt);
                z        = IP*rho;
            end    
        end   
   end
    
end