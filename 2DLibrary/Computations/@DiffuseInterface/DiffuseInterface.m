classdef DiffuseInterface < Computation
    
   properties (Access = public)              
       IC                     
       IntSubArea
              
       phi = [],uv  = []
       theta=[]
       errors = []
       
       IsolineInterfaceY2=[],StagnationPoint=[]
       %surfaceTension = 4/3;
       
       filename
   end
   
   methods (Abstract = true, Access = public)
       [A,b] = FullStressTensorIJ(this,phi,i,j)
       [A,b] = ContinuityMomentumEqs(this,phi)
       [bulkError,bulkAverageError] = DisplayFullError(this) 
   end
   
   methods (Access = public)          
        function this = DiffuseInterface(config)  
            this@Computation(config);            
            if(~isfield(this.optsNum.PhysArea,'y1Max'))
                this.optsNum.PhysArea.y1Max = inf;
            end      
            if(~isfield(this.optsNum.PhysArea,'NBorder'))
                this.optsNum.PhysArea.NBorder = [];
            end
            
            if(~isfield(this.optsNum.PlotArea,'y2Min'))
                this.optsNum.PlotArea.y2Min = this.optsNum.PhysArea.y2Min;
            end
            if(~isfield(this.optsNum.PlotArea,'y2Max'))
                this.optsNum.PlotArea.y2Max = this.optsNum.PhysArea.y2Max;
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

            bxArea          = struct('y1Min',this.optsNum.PhysArea.IntInterval(1),...
                                     'y1Max',this.optsNum.PhysArea.IntInterval(2),...
                                     'y2Min',this.optsNum.PhysArea.y2Min,...
                                     'y2Max',this.optsNum.PhysArea.y2Max,...
                                     'N',[100,100]);
            
            BX              = Box(bxArea);
            IntBx           = BX.ComputeIntegrationVector();
            this.IntSubArea = IntBx*this.IC.SubShapePts(BX.GetCartPts());
                     
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
        function spt = FindStagnationPoint(this,iguess1,iguess2)
            
            uv = this.uv;
            
            if(nargin < 2)
                iguess1 = [-2,2];
            end
            if(nargin < 3)
                iguess2 = [-2,this.IC.y2Max-2];
            end
            
            fsolveOpts = optimset('Display','off');
            [spt1,~,flag1] = fsolve(@ValueAtXY,iguess1,fsolveOpts);                        
            if(flag1 >= 1)
                spt.y1_kv = spt1(1);
                spt.y2_kv = spt1(2);
                disp(['Stagnation point found at (',num2str(spt1(1)),',',num2str(spt1(2)),')']);
            end
                        
            [spt2,~,flag2] = fsolve(@ValueAtXY,iguess2,fsolveOpts);
            if(flag2 >= 1)
                spt.y1_kv = [spt.y1_kv;spt2(1)];
                spt.y2_kv = [spt.y2_kv;spt2(2)];
                disp(['Stagnation point found at (',num2str(spt2(1)),',',num2str(spt2(2)),')']);
            end
            
            if((flag1 < 1) && (flag2 < 1))
                disp('No stagnation point found');
                spt = [];
            else
                this.StagnationPoint = spt;                
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
        function PlotResultsMu(this)                         
            figure('Position',[0 0 800 600],'color','white');
            this.IC.doPlots(GetMu(this),'contour');             
            PlotU(this); hold on;             
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
            if(~isempty(this.theta))
                this.PlotSeppecherSolution(this.theta,this.phi);
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
        function AddStreamlines(this)
            for i = 1:3
                [y10,y20] = ginput(1);   
                this.IC.doPlotsStreamlines(this.uv,y10,y20); %IC.doPlotsFlux(u_flow)(mu);
            end
        end
        function PlotU(this,uv,y1Pts,y2Pts) 
            
            n = 20;
            if((nargin<2) || isempty(uv))
                if(isempty(this.uv))
                    return;
                else
                    uv = this.uv;
                end
            end
            y2Min = this.optsNum.PhysArea.y2Min;
            y2Max = this.optsNum.PhysArea.y2Max;
            
            y1Min = this.optsNum.PlotArea.y1Min;
            y1Max = this.optsNum.PlotArea.y1Max;
                        
            y2L = y2Min + (y2Max-y2Min)*(0:n-1)'/(n-1);            
            y1L = y1Min + (y1Max-y1Min)*(0:n-1)'/(n-1);            
            
            startPtsy1    = [y1Max*ones(size(y2L))-0.1;...
                         y1Min*ones(size(y2L))+0.1;...
                         y1L];
            startPtsy2    = [y2L;y2L;y2Max*ones(size(y1L))];
            
            if(nargin >= 4)
               startPtsy1 = [startPtsy1;y1Pts]; 
               startPtsy2 = [startPtsy2;y2Pts];
            end
            
            this.IC.doPlotsStreamlines(uv,startPtsy1,startPtsy2); %IC.doPlotsFlux(u_flow)(mu);
            hold on;
            this.IC.doPlotsFlux(uv);
            
            sp = this.StagnationPoint;
            if(~isempty(sp))          
                hold on;
                plot(sp.y1_kv,sp.y2_kv,'or','MarkerFaceColor','r','MarkerSize',10); 
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
        function SavePlotResults(this)                                
            
            PlotResultsPhi(this);
            SaveCurrentFigure(this,[this.filename '_Density']);            
            
            PlotResultsMu(this);
            SaveCurrentFigure(this,[this.filename '_ChemPot']);            
        end        
        function PlotErrorIterations(this)           
            disp('*** Results ***');
            if(~isempty(this.theta))
                disp(['theta = ',num2str(this.theta*180/pi),' [deg]']);
            end
            if(isfield(this.errors,'aIter'))
                disp(['a = ',num2str(this.errors.aIter(end))]);
            end            

            legStr = {};
            figure('color','white');
            if(isfield(this.errors,'errorIterations'))
                disp(['Max error of equations excluding boundaries: ',num2str(this.errors.errorIterations(end))]);
                semilogy(this.errors.errorIterations,'ro','MarkerFaceColor','r'); hold on;
                legStr{end+1} = 'Error';
            end
            if(isfield(this.errors,'errorAverage'))
                semilogy(this.errors.errorAverage,'ko','MarkerFaceColor','k');
                legStr{end+1} = 'Average error';
            end
            if(length(legStr) > 1)
                legend(legStr);            
            end
            xlabel('Iteration');
            ylabel('Error');
                        
            SaveCurrentFigure(this,[this.filename '_ErrorIterations']);                 
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
        function theta  = FindInterfaceAngle(this,phi)

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
        
        function AnalyzeScalarQuantity(this,f,interval)
            if(nargin == 2)
                interval = [-10,10];
            end
            y2Max = this.optsNum.PhysArea.y2Max;
            noCuts = 5;
            optsC  = {'b','g','m','k','r'};
            
            figure('Position',[0 0 1000 800]);
            subplot(2,1,1);
            leg = {};
            for i= 0:1:(noCuts-1)
                y2 = y2Max*i/(noCuts-1);
                this.IC.doPlotFLine(interval,[y2 y2],f,[],optsC{i+1});  hold on;
                
                IP       = this.IC.SubShapePtsCart(struct('y1_kv',-inf,'y2_kv',y2));
                plot(interval,[1,1]*(IP*f),[optsC{i+1},'--']); hold on;
                IP       = this.IC.SubShapePtsCart(struct('y1_kv',inf,'y2_kv',y2));
                plot(interval,[1,1]*(IP*f),[optsC{i+1},'-.']);  hold on;
                
                leg{end+1} = ['y2 = ',num2str(y2)];
                leg{end+1} = '';
                leg{end+1} = ['y2 = ',num2str(y2),'y1=-inf'];
                leg{end+1} = ['y2 = ',num2str(y2),'y1=inf'];
            end
            legend(leg,'Location','eastoutside');
            
            subplot(2,1,2);
            this.IC.doPlotFLine([-inf,-inf],[0 y2Max],f,[],'r'); hold on;
            this.IC.doPlotFLine([inf,inf],[0 y2Max],f,[],'b'); hold on;
            
            this.IC.doPlotFLine([-4 -4],[0 y2Max],f);  hold on;
            this.IC.doPlotFLine([0 0],[0 y2Max],f);  hold on;
            this.IC.doPlotFLine([4 4],[0 y2Max],f);  hold on;
            
            legend({'-inf','inf'},'Location','eastoutside');
        end
        
        [phi,muDelta]  = GetEquilibriumDensityR(this,mu,theta,nParticles,phi,ptC)                  
        [phi,uv]       = SolveFull(this,ic)
        
        
        function [v_mu,A_mu] = ChemicalPotential(this,uv,phi,G)
            %[uv;phi;G]     
            Cn             = this.optsPhys.Cn;
            Ind            = this.IC.Ind;
            
            M              = this.IC.M;
            F              = false(M,1);   
            T              = true(M,1);            
            Z              = zeros(M);
            
            Diff           = this.IC.Diff;            
            
    
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            % **************
            
            A_mu           = [Z,Z,...
                              diag(fWPP)-Cn*Diff.Lap,...
                              -eye(M)];
            v_mu           = fWP - Cn*Diff.Lap*phi - G;
            
            % (BC4.a) nu*grad(phi) = 0            
            A_mu(Ind.bottom|Ind.top,:)         = 0;
            A_mu(Ind.bottom|Ind.top,[F;F;T;F]) = Diff.Dy2(Ind.bottom|Ind.top,:);    
            v_mu(Ind.bottom|Ind.top)           = Diff.Dy2(Ind.bottom|Ind.top,:)*phi;
       
        end
   end
    
end