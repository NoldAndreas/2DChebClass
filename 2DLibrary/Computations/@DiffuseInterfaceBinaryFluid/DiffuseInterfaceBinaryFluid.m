classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p  % pressure       
       mu % chemical potential       
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

             G    = this.mu;  
             vec  = [this.uv;this.phi;G;this.p];
       end
       
       function SavePlotResults(this)           
           SavePlotResults@DiffuseInterface(this);
           
           figure('Position',[0 0 800 600],'color','white');
           PlotU(this); hold on;            
           this.IC.doPlots(this.p,'contour');     
           SaveCurrentFigure(this,[this.filename '_Pressure']);                
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
       
       
       function [v_mu,A_mu] = ChemicalPotential(this,uv,phi,G)
           [v_mu,A_mu] = ChemicalPotential@DiffuseInterface(this,uv,phi,G);
           A_mu = [A_mu,zeros(this.IC.M)];
       end
       function [v_G,A_G] = Phasefield(this,uv,phi,G)
           
            m              = this.optsPhys.mobility;
            Ind            = this.IC.Ind;
            Diff           = this.IC.Diff;
            M    = this.IC.M;
            EYM  = eye(M); 
            Z    = zeros(M);
            F       = false(M,1);   T       = true(M,1);    
            
            uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];        
            gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];
   
           
            A_G            = [-gradPhiMT,-uvTdiag*Diff.grad,m*Diff.Lap,Z];
            v_G            = m*Diff.Lap*G - gradPhiMT*uv;
           
            % (BC4) nu*grad(G) = 0        
            A_G(Ind.top|Ind.bottom,:)           = 0;
            A_G(Ind.top|Ind.bottom,[F;F;F;T;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
            v_G(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*G;                                    


            %(BC7) 
            indBC7 = Ind.left | Ind.right;
            A_G(indBC7,:)           = 0;
            A_G(indBC7,[F;F;F;T;F]) = Diff.Dy2(indBC7,:);
            v_G(indBC7,:)           = Diff.Dy2(indBC7,:)*G;

            %(BC8) [uv;phi;G;p]        
            hM  = - m*this.IC.borderBottom.IntSc*Diff.DDy2;        
            A_G(Ind.right & Ind.bottom,:)           = 0;
            A_G(Ind.right & Ind.bottom,[F;T;F;F;F]) = this.IC.borderBottom.IntSc*(diag(phi)*Diff.Dy2);
            A_G(Ind.right & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) ...
                                                      - EYM(Ind.left & Ind.bottom,:) ...                                              
                                                      + this.IC.borderBottom.IntSc*(Diff.Dy2*(uv(1+end/2:end)));
            A_G(Ind.right & Ind.bottom,[F;F;F;T;F]) = hM;
            v_G(Ind.right & Ind.bottom)             = phi(Ind.right & Ind.bottom) ...
                                                      - phi(Ind.left & Ind.bottom) ...
                                                      + this.IC.borderBottom.IntSc*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                                      + hM*G;
       end
       
   end
end