classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p  % pressure       
       mu % chemical potential
       s = 0.
   end
       
	methods (Access = public) 
       function this = DiffuseInterfaceBinaryFluid(config)           
             this@DiffuseInterface(config);
       end
       
       IterationStepFullProblem(this,vec)
       [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi,findTheta)
       [rho,theta]         = GetEquilibriumDensity_Flux(this,theta,phi,uv)
       
       [A,b] = ContinuityMomentumEqs(this,phi,mu)
       [A,b] = ContinuityMomentumEqs_mu_p_uv(this,phi)
       
       SolveMovingContactLine(this,noIterations)
       [p,uv,A,b,a] = GetVelocityAndChemPot(this,phi,mu,theta)
       [p,mu,uv,A,b,a] = GetVelocityPressureAndChemPot(this,phi,theta)
       
       mu = SolvePhasefieldForChemPot(this,uv,phi)
       
       function [bulkError,bulkAverageError] = DisplayFullError(this)
            M        = this.IC.M;
            m        = this.optsPhys.mobility;
            Cn       = this.optsPhys.Cn;
            Diff     = this.IC.Diff;
            
            uv       = this.uv;
            phi      = this.phi;
            p        = this.p;
            mu       = this.mu;
            
            [Af,bf]  = ContinuityMomentumEqs(this,phi,mu);                                     
            Convect    = [diag(uv(1:end/2)),diag(uv(1+end/2:end))]*this.IC.Diff.grad/m;
            
            error    = Af*[p;uv] - bf;
            
            errorPhi    =  - Convect*phi + Diff.Lap*mu;
                    %absGradRho  = ( (Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2 );
                            %+ 12*Cn*phi.*absGradRho ...
                            %+ 2*Cn*(3*phi.^2 -1).*(Diff.Lap*phi)...
                            %- Cn*Diff.Lap2*phi;            
            
            PrintErrorPos(error(1:M),'continuity equation',this.IC.Pts);            
            PrintErrorPos(error(1+M:2*M),'y1-momentum equation',this.IC.Pts);            
            PrintErrorPos(error(2*M+1:end),'y2-momentum equation',this.IC.Pts);                        
            PrintErrorPos(errorPhi,'Phasefield equation',this.IC.Pts);
            
            PrintErrorPos(error(repmat(~this.IC.Ind.bound,3,1)),'bulk continuity and momentum equations');               
            PrintErrorPos(errorPhi(~this.IC.Ind.bound),'bulk phasefield equations');               
            
            bulkError        = max(abs(error(repmat(~this.IC.Ind.bound,2,1))));
            bulkAverageError = mean(abs(error(repmat(~this.IC.Ind.bound,2,1))));
       end        
       function vec = GetInitialCondition(this) 
           
            if(~isempty(this.uv) && ...
                ~isempty(this.p) && ...
                ~isempty(this.mu) && ...
                ~isempty(this.phi))
                vec  = [this.uv;this.phi;this.mu - this.s*this.phi;this.p];
                return;
            end
           
            M = this.IC.M;
            if(isempty(this.theta))
                otherInput.thetaInitialGuess = pi/2;
            else
                otherInput.thetaInitialGuess = this.theta;
            end            

             if(otherInput.thetaInitialGuess == pi/2)
                theta  = pi/2;
                phi    = InitialGuessRho(this);
                mu     = 0;
             else
                theta  = otherInput.thetaInitialGuess;
                phi    = this.phi;
                mu     = this.GetMu();
             end

             this.uv  = GetBoundaryCondition(this);%,theta,phi);            
             this.p   = zeros(M,1);
             this.mu  = zeros(M,1);
             this.phi = phi;                   

             G    = this.mu - this.s*this.phi;             
             vec  = [this.uv;this.phi;G;this.p];
       end
       
       function SavePlotResults(this)
           global dirData
           SavePlotResults@DiffuseInterface(this);
           
           figure('Position',[0 0 800 600],'color','white');
           PlotU(this); hold on;            
           this.IC.doPlots(this.p,'contour');     
           print2eps([dirData filesep this.filename '_Pressure'],gcf);
           saveas(gcf,[dirData filesep this.filename '_Pressure.fig']);
       end
       
        function PlotU(this) 
            
            if(isempty(this.StagnationPoint))
               FindStagnationPoint(this,[-2,2],[-2,this.IC.y2Max-2]); 
            end            
            
            y1Pts = [0;...
                     this.optsNum.PlotArea.y1Min+0.1;...
                     this.optsNum.PlotArea.y1Max-0.1];
            y2Pts = [this.optsNum.PlotArea.y2Max/2;...
                     0.5;...
                     this.optsNum.PlotArea.y2Max-0.5];
            
            PlotU@DiffuseInterface(this,[],y1Pts,y2Pts);
        end

       
       function CheckResultResolution(this)
           figure('Position',[0 0 1000 1000],'color','white');
           
           y2M = this.IC.y2Max;
           
           subplot(2,2,1);           
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.mu);           
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.mu);
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.mu);
           title('$\Delta \mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,2);
           this.IC.doPlotFLine([-20 20],[0 0],this.mu);           
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.mu);
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.mu);
           title('$\mu$','Interpreter','Latex','fontsize',15);
           
           subplot(2,2,3);
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.uv(1:end/2));
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.uv(1:end/2));
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.uv(1:end/2));
           title('$\Delta u$','Interpreter','Latex','fontsize',15);           
           
           subplot(2,2,4);
           this.IC.doPlotFLine([-20 20],[0 0],this.IC.Diff.Lap*this.uv(1+end/2:end));
           this.IC.doPlotFLine([-20 20],[y2M y2M]/2,this.IC.Diff.Lap*this.uv(1+end/2:end));
           this.IC.doPlotFLine([-20 20],[y2M y2M],this.IC.Diff.Lap*this.uv(1+end/2:end));      
           title('$\Delta v$','Interpreter','Latex','fontsize',15);
           
       end
       
   end
end