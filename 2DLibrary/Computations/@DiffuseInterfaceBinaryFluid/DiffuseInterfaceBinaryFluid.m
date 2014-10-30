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
                vec  = [this.uv;this.phi;this.mu;this.p];
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
       
       
       function [v_cont,A_cont] = Continuity(this,uv,phi,G,p)            
           
           Ind  = this.IC.Ind;
           Diff = this.IC.Diff;
           M    = this.IC.M;
           EYM  = eye(M);
           Z    = zeros(M);
           
           F    = false(M,1);   
           T    = true(M,1);    
           
           Cn          = this.optsPhys.Cn;
           Cak         = this.optsPhys.Cak;
           y2Max       = this.optsNum.PhysArea.y2Max;
           
           phiM        = phi(Ind.left & Ind.top);
           phiP        = phi(Ind.right & Ind.top);  
        
           pM          = p(Ind.left & Ind.top);
           pP          = p(Ind.right & Ind.top);
        
           IntPathUpLow   = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc;     
           
           A_cont         = [Diff.div,Z,Z,Z];    
           v_cont         = Diff.div*uv;
            
           % (BC1) p = 0  at y1 = +/- infinity
           A_cont(Ind.left,:)           = 0;
           A_cont(Ind.left,[F;F;F;F;T]) = EYM(Ind.left,:);
           v_cont(Ind.left)             = p(Ind.left);

            % (BC2) [uv;phi;G;p]

           [fpW_mInf,fW_mInf] = DoublewellPotential(phiM,Cn);
           [fpW_pInf,fW_pInf] = DoublewellPotential(phiP,Cn);

           indBC2 = Ind.right & Ind.top;

           A_cont(Ind.right,:)           = 0;
           A_cont(Ind.right,[F;F;F;F;T]) = Diff.Dy2(Ind.right,:);
           v_cont(Ind.right,:)           = Diff.Dy2(Ind.right,:)*p;

           A_cont(indBC2,[T;T;F;F;F]) = IntPathUpLow*[Diff.Dy2 , Diff.Dy1];              

           A_cont(indBC2,[F;F;T;F;F])                 = -Cn/Cak*IntPathUpLow*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1);
           A_cont(indBC2,[F;F;Ind.top&Ind.right;F;F]) = ...
                            A_cont(indBC2,[F;F;Ind.top&Ind.right;F;F]) + ...
                            y2Max*fpW_pInf/Cak;
           A_cont(indBC2,[F;F;Ind.top&Ind.left;F;F])  = ...
                            A_cont(indBC2,[F;F;Ind.top&Ind.left;F;F]) - ...
                            y2Max*fpW_mInf/Cak;

           A_cont(indBC2,[F;F;F;F;Ind.top&Ind.right]) = -y2Max;  
           A_cont(indBC2,[F;F;F;F;Ind.top&Ind.left])  = y2Max;  


           v_cont(indBC2) = IntPathUpLow*([Diff.Dy2 , Diff.Dy1]*uv ...
                                  - Cn/Cak*((Diff.Dy1*phi).*(Diff.Dy2*phi)))...
                            +y2Max*((-pP + fW_pInf/Cak) - (-pM + fW_mInf/Cak));

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
            A_G(Ind.left & Ind.bottom,:)           = 0;
            A_G(Ind.left & Ind.bottom,[F;T;F;F;F]) = this.IC.borderBottom.IntSc*(diag(phi)*Diff.Dy2);
            A_G(Ind.left & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) ...
                                                      - EYM(Ind.left & Ind.bottom,:) ...                                              
                                                      + this.IC.borderBottom.IntSc*(Diff.Dy2*(uv(1+end/2:end)));
            A_G(Ind.left & Ind.bottom,[F;F;F;T;F]) = - m*this.IC.borderBottom.IntSc*Diff.DDy2;
            v_G(Ind.left & Ind.bottom)             = phi(Ind.right & Ind.bottom) ...
                                                      - phi(Ind.left & Ind.bottom) ...
                                                      + this.IC.borderBottom.IntSc*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                                      - m*this.IC.borderBottom.IntSc*Diff.DDy2*G;
                              
            % ****************
            IP                      = this.IC.SubShapePtsCart(this.RightCapillary.GetCartPts);           
            IP0                     = this.IC.SubShapePtsCart(struct('y1_kv',0,'y2_kv',0));
            IntPath_Half_Low        = this.RightCapillary.borderBottom.IntSc*IP;            
            
            A_G(Ind.right & Ind.bottom,:)           = 0;
            A_G(Ind.right & Ind.bottom,[F;T;F;F;F]) = IntPath_Half_Low*(diag(phi)*Diff.Dy2);
            A_G(Ind.right & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) - IP0 ...                                              
                                                      + IntPath_Half_Low*(Diff.Dy2*(uv(1+end/2:end)));
            A_G(Ind.right & Ind.bottom,[F;F;F;T;F]) = - m*IntPath_Half_Low*Diff.DDy2;
            
            v_G(Ind.right & Ind.bottom) = phi(Ind.right & Ind.bottom) - IP0*phi ...
                                          + IntPath_Half_Low*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                          - m*IntPath_Half_Low*Diff.DDy2*G;

       end
       function [v_mom,A_mom]   = Momentum(this,uv,phi,G,p)                
           
            %[uv;phi;G;p]         
            Diff = this.IC.Diff;
            Cak  = this.optsPhys.Cak;            
           
            A_mom          = [Diff.LapVec,...
                               diag(repmat(G,2,1))*Diff.grad/Cak,...                               
                               [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)]/Cak,...
                               - Diff.grad];
            v_mom          = Diff.LapVec*uv ...
                              + repmat(G,2,1).*(Diff.grad*phi)/Cak ...
                              - Diff.grad*p;
       end
   end
end