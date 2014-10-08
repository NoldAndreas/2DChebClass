classdef DiffuseInterface < handle
    
   properties (Access = public)              
       IC              
       optsNum,optsPhys    
       IntSubArea
       PlotBCShape
       
       phi = [],uv  = [],theta=[],a
       errors = []
       
       IsolineInterfaceY2=[],StagnationPoint=[]
       %surfaceTension = 4/3;
       
       configName,filename
   end
   
   methods (Abstract = true, Access = public)
       [A,b] = FullStressTensorIJ(this,phi,i,j)
       [A,b] = ContinuityMomentumEqs(this,phi)
       [bulkError,bulkAverageError] = DisplayFullError(this) 
   end
   
   methods (Access = public)          
        function this = DiffuseInterface(config)             
            this.optsNum         = config.optsNum;
            this.optsPhys        = config.optsPhys;                                                               
            this.configName    = SaveConfig(config,'Configurations');                                                                         
            if(~isfield(this.optsNum.PhysArea,'y1Max'))
                this.optsNum.PhysArea.y1Max = inf;
            end            
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
                   'gradLap' ;'gradDiv'; 'LapVec';'LapDiv';'Lap2'};
            this.IC.ComputeDifferentiationMatrix(Sel);
            
            this.IC.SetUpBorders(this.optsNum.PhysArea.NBorder);            
            
            this.IC.Ind.fluidInterface = [];
            if(isfield(this.optsPhys,'fluidInterface'))
                this.IC.Ind.fluidInterface = this.IC.Ind.(this.optsPhys.fluidInterface);
            end

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
        function phi = InitialGuessRho(this)
            PtsCart    = this.IC.GetCartPts();
            Cn         = this.optsPhys.Cn;
            
            phi        = tanh(PtsCart.y1_kv/Cn);
        end                                                  
                
        deltaX                      = GetDeltaX(this,phi,theta)              
        function [a,deltaX]         = Get_a_deltaX(this,phi,theta)
            
            Ind              = this.IC.Ind;
            UWall            = this.optsPhys.UWall;    
            InterpOntoBorder = this.IC.borderTop.InterpOntoBorder;
            ptsBorderTop     = this.IC.borderTop.Pts;
            y2Max            = this.optsNum.PhysArea.y2Max;
            phiL             = phi(Ind.top & Ind.left);
            phiR             = phi(Ind.top & Ind.right);

            
            IntNormal_Path   = this.IC.borderTop.IntNormal_Path;       
            phi_m            = this.optsPhys.phi_m;    
            
            deltaX              = GetDeltaX(this,phi,theta);
            ptsBorderTop.y1_kv  = ptsBorderTop.y1_kv - deltaX;                        
            u_flow              = GetSeppecherSolutionCart(ptsBorderTop,UWall,0,0,theta);                                    

            phiBorder  = InterpOntoBorder*phi;
            phiBorder2 = repmat(phiBorder,2,1);

            massFlux   = ((phiR+phi_m)-(phiL+phi_m))*(y2Max-0)*UWall;      %due to mapping to infinity

            fsolveOpts = optimset('Display','off');
            [a,~,flag] = fsolve(@massInflux,0,fsolveOpts);            
            if(flag < 1)
                cprintf('*r','CorrectVelocityProfile: Finding parameter a: No solution found.\n');                                
            else
                disp(['a = ',num2str(a)]);
            end         
            
            function m = massInflux(a)
                u_corrected = u_flow .*(1 + a*(phiL-phiBorder2).^2.*(phiR-phiBorder2).^2);
                m          = IntNormal_Path*(u_corrected.*(phiBorder2+phi_m)) + massFlux;    
            end       
        end                    
        function [uvBound,a,deltaX] = GetBoundaryCondition(this,theta,phi)
            UWall            = this.optsPhys.UWall;    
            PtsCart          = this.IC.GetCartPts();
            Ind              = this.IC.Ind;            
            y2_kv            = PtsCart.y2_kv;
            y2Min            = this.optsNum.PhysArea.y2Min;
            y2Max            = this.optsNum.PhysArea.y2Max;
            M                = this.IC.M;
            
            if(length(UWall) == 1)
                phiL             = phi(Ind.top & Ind.left);
                phiR             = phi(Ind.top & Ind.right);
                
                [a,deltaX] = Get_a_deltaX(this,phi,theta);    
                u_flow     = GetSeppecherSolutionCart([PtsCart.y1_kv - deltaX,...
                                         PtsCart.y2_kv],UWall,0,0,theta);          
                phiBorder2 = repmat(phi,2,1);
                uvBound    = u_flow .*(1 + a*(phiL-phiBorder2).^2.*(phiR-phiBorder2).^2);   
            else
                a = [];  deltaX = [];
                uvBound    = [UWall(1) + ...
                            (UWall(2)-UWall(1))*(y2_kv-y2Min)/(y2Max-y2Min);...
                            zeros(M,1)];                
            end
        end
                         
                
        function spt = FindStagnationPoint(this)
            
            uv = this.uv;
            
            fsolveOpts = optimset('Display','off');
            [spt,~,flag] = fsolve(@ValueAtXY,[-2,2],fsolveOpts);                        
            
            if(flag < 1)
                disp('No stagnation point found');
                spt = [];
            else
                this.StagnationPoint = spt;
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
        [A,b] = Div_FullStressTensor(this,phi)
        function mu_s = GetMu(this,phi)
            if(nargin == 1)
                phi = this.phi;
            end
            Cn       = this.optsPhys.Cn;
            mu_s     = DoublewellPotential(phi,Cn) - Cn*(this.IC.Diff.Lap*phi);               
        end                   
        
        %Analysis functions
        function ComputeInterfaceContour(this)
            y2           = this.IC.Pts.y2;
            phi          = this.phi; 

            fsolveOpts   = optimset('Display','off');
            interface    = zeros(size(y2));

            y1Bottom = this.IC.GetCartPts.y1_kv(this.IC.Ind.bottom);
            [~,j]    = min(abs(phi(this.IC.Ind.bottom)));
            y1I      = y1Bottom(j);
                        
            for i = 1:length(y2)
                pt.y2_kv     = y2(i);        
                [interface(i),~,flag] = fsolve(@phiX1,y1I,fsolveOpts);        
                if(flag < 1)
                    cprintf('*r',['ComputeInterfaceContour: Interface not found at y2 = ',num2str(y2(i)),'\n']);
                    return;
                end
                y1I          = interface(i);        
            end

            this.IsolineInterfaceY2 = interface;    

            function z = phiX1(y1)
                pt.y1_kv = y1;
                IP = this.IC.SubShapePtsCart(pt);
                z  = IP*phi;
            end   
        end
        
        %Plotting                           
        function PlotResultsMu(this,mu,uv) 
            if(nargin == 1)
                uv  = this.uv;
                mu  = GetMu(this,this.phi);            
            end
            
            figure('Position',[0 0 800 600],'color','white');
            this.IC.doPlots(mu,'contour');             
            PlotU(this,uv); hold on;             
        end
        function PlotResultsPhi(this)
            figure('Position',[0 0 800 600],'color','white');
            
            PlotU(this); hold on;             
            this.IC.doPlots(this.phi,'contour');     
                        
            hold on;
            if(~isempty(this.IsolineInterfaceY2))
                plot(this.IsolineInterfaceY2,this.IC.Pts.y2,...
                                                    'k','linewidth',3);
            end
            if(sum(this.IC.Ind.fluidInterface)>0)
                this.PlotSeppecherSolution(theta,phi);
            end
        end        
        function PlotSeppecherSolution(this,theta,phi)
            if(nargin == 1)                
                phi = this.phi;            
                theta = this.theta;
            end
            UWall   = this.optsPhys.UWall;
%            D_A     = this.optsPhys.D_A;                        
%             u_flow = GetSeppecherSolutionCart(this.PlotBCShape.GetCartPts,...
%                                               UWall,D_A,D_B,theta);                                                                        
%             this.PlotBCShape.doPlots(u_flow,'flux',struct('reshape',false,'linecolor','m'));
            
            PtsCart = this.IC.GetCartPts(); 
            if(nargin > 2)
                deltaX = GetDeltaX(this,phi,theta);             
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
%                 subplot(1,2,2);  this.IC.doPlotFLine([-100,100],[y2Max,y2Max],phi.*u_flow(end/2+1:end),'CART'); 
%             else
%                 figure('Position',[0 0 800 800],'Color','white');            
%                 PlotU(this,u_flow);            
%             end
%             title('Check accuracy of map');
        end          
        function SavePlotResults(this)
            global dirData                        
            
            PlotResultsPhi(this);
            print2eps([dirData filesep this.filename '_Density'],gcf);
            saveas(gcf,[dirData filesep this.filename '_Density.fig']);
            
            PlotResultsMu(this);
            print2eps([dirData filesep this.filename '_ChemPot' ],gcf);
            saveas(gcf,[dirData filesep this.filename '_ChemPot.fig']);
        end        
        function PlotU(this,uv) 
            
            if(nargin<2)
                if(isempty(this.uv))
                    return;
                else
                    uv = this.uv;
                end
            end
            y2Max = this.optsNum.PhysArea.y2Max;
            
            y1Min = this.optsNum.PlotArea.y1Min;
            y1Max = this.optsNum.PlotArea.y1Max;
            
            y2L = (1:y2Max)';
            y1L = (y1Min:y1Max)';            
            
            startPtsy1    = [y1Max*ones(size(y2L))-0.1;...
                         y1Min*ones(size(y2L))+0.1;...
                         y1L];
            startPtsy2    = [y2L;y2L;y2Max*ones(size(y1L))];
            this.IC.doPlotsStreamlines(uv,startPtsy1,startPtsy2); %IC.doPlotsFlux(u_flow)(mu);
            
            sp = this.StagnationPoint;
            if(~isempty(sp))          
                hold on;
                plot(sp(1),sp(2),'or','MarkerFaceColor','r','MarkerSize',10); 
                hold on;
            end
        end           
        function PlotInterfaceAnalysis(this)            
            
            y2      = this.IC.Pts.y2;
            L       = this.optsNum.PhysArea.y2Max;
            
            Diff    = barychebdiff(y2,1);
            thetaY2 = pi/2 + atan(Diff.Dx*this.IsolineInterfaceY2);%*180/pi;
            
            A = -0.34; k = 0.2; 
            analytic_BriantYeomans = A* (-(1+exp(-k*L)) + (exp(-k*y2) + exp(-k*(L-y2))));            
            
            %maxthetaY2 = max(ab
            
            figure('color','white','Position',[0 0 800 800]);
            
            plot(y2,cos(thetaY2),'k','linewidth',2); hold on;
            plot(y2,analytic_BriantYeomans,'k--','linewidth',2);
            
            xlabel('$y_2$','Interpreter','Latex','fontsize',20);
            ylabel('$cos(\theta)$','Interpreter','Latex','fontsize',20);
        end

        function PlotErrorIterations(this)           
            disp('*** Results ***');
            disp(['theta = ',num2str(this.theta*180/pi),' [deg]']);
            disp(['a = ',num2str(this.errors.aIter(end))]);
            disp(['Max error of equations excluding boundaries: ',num2str(this.errors.errorIterations(end))]);

            figure('color','white');
            semilogy(this.errors.errorIterations,'ro','MarkerFaceColor','r'); hold on;
            semilogy(this.errors.errorAverage,'ko','MarkerFaceColor','k');
            legend({'Maximal error','Average error'});
            
            xlabel('Iteration');
            ylabel('Error');
            
            global dirData
            
            print2eps([dirData filesep this.filename '_ErrorIterations'],gcf);
            saveas(gcf,[dirData filesep  this.filename '_ErrorIterations.fig']);
        end
        
        %Old
        D_B         = SetD_B(this,theta,phi,initialGuessDB)
        function Interp = ResetOrigin(this,phi)
            
            pt.y2_kv  =  0;
            DeltaY1   = fsolve(@phiX1,-1);%,fsolveOpts);
            
            ptsCartShift       = this.IC.GetCartPts();
            ptsCartShift.y1_kv = ptsCartShift.y1_kv + DeltaY1;
            
            Interp = this.IC.SubShapePtsCart(ptsCartShift);                       
            
            function z = phiX1(y1)
                pt.y1_kv = y1;
                IP       = this.IC.SubShapePtsCart(pt);
                z        = IP*phi;
            end    
        end
        function theta = FindInterfaceAngle(this,phi)

            y2M = 0; y2P = 10;
          %  fsolveOpts   = optimset('Display','off');                

            pt.y2_kv  =  y2M;
            y1CartStart = fsolve(@phiX1,-1);%,fsolveOpts);

            pt.y2_kv  = y2P;
            y1CartEnd = fsolve(@phiX1,2);%,fsolveOpts);                

            alpha = atan((y1CartStart-y1CartEnd)/(y2P-  y2M));
            theta = alpha + pi/2;
            
            disp(['Contact angle: ',num2str(theta*180/pi),'[deg].']);

            function z = phiX1(y1)
                pt.y1_kv = y1;
                IP       = this.IC.SubShapePtsCart(pt);
                z        = IP*phi;
            end    
        end   
        
        [phi,muDelta] = GetEquilibriumDensityR(this,mu,theta,nParticles,phi,ptC)                  
        [phi,uv]       = SolveFull(this,ic)
   end
    
end