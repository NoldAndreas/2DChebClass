classdef DiffuseInterfaceSingleFluid < DiffuseInterface
	properties (Access = public)              
        mu
        a=[],deltaX=[]
    end
         
    methods (Access = public) 
         function this = DiffuseInterfaceSingleFluid(config)           
             this@DiffuseInterface(config);
         end
                 
         SolveMovingContactLine(this,maxIterations)       
         [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi,findTheta)
         [mu,uv,A,b,a]       = GetVelocityAndChemPot(this,phi,theta)
         [A,b]               = ContinuityMomentumEqs(this,phi)
                        
         [A,b]       = FullStressTensorIJ(this,phi,i,j)           
         function p = GetPressure_from_ChemPotential(this,mu,phi_ig)
            Cn     = this.optsPhys.Cn;
            phi_m  = this.optsPhys.phi_m;
            %for bulk
            % 0 = W' - m
            % p = - W + mu*(phi + phi_m)
            
            fsolveOpts = optimset('Display','off');
            [phi,~,exitflag]   = fsolve(@f,phi_ig,fsolveOpts);
            if(exitflag < 1)
                cprintf('*r','No solution for density for given chemical potential found');
            end
            
            [~,W] = DoublewellPotential(phi,Cn);
            p     = - W + mu*(phi+phi_m);
            
            function y = f(phi)
                y = DoublewellPotential(phi,Cn) - mu;
            end
            
         end  
        
         IterationStepFullProblem(this,noIterations)
         IterationStepFullProblem_Seppecher(this,noIterations)
         function vec = GetInitialCondition(this) 
           
            if(~isempty(this.uv) && ...                
                ~isempty(this.mu) && ...
                ~isempty(this.phi))
                vec  = [this.uv;this.phi;this.mu];
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

             this.uv  = GetBoundaryCondition(this,theta,phi);                         
             this.mu  = zeros(M,1);
             this.phi = phi;
             
             vec  = [this.uv;this.phi;this.mu];
         end
        
         function [bulkError,bulkAverageError] = DisplayFullError(this,phi,uv)            
            M        = this.IC.M;
            
            [Af,bf]  = ContinuityMomentumEqs(this,phi);                         
            mu_s     = GetMu(this,phi);
            
            error    = Af*[mu_s;uv] - bf;
            PrintErrorPos(error(1:M),'continuity equation',this.IC.Pts);            
            PrintErrorPos(error(1+M:2*M),'y1-momentum equation',this.IC.Pts);            
            PrintErrorPos(error(2*M+1:end),'y2-momentum equation',this.IC.Pts);                        
            
            PrintErrorPos(error(repmat(~this.IC.Ind.bound,2,1)),'bulk continuity and momentum equations');               
            
            bulkError        = max(abs(error(repmat(~this.IC.Ind.bound,2,1))));
            bulkAverageError = mean(abs(error(repmat(~this.IC.Ind.bound,2,1))));
         end        
         function [v_cont,A_cont] = Continuity(this,uv,phi,G)
            
            Cn             = this.optsPhys.Cn;
            Cak            = this.optsPhys.Cak; 
            Ind            = this.IC.Ind;
            zeta           = this.optsPhys.zeta;
                        
            lbC            = Ind.left & Ind.bottom;
            rbC            = Ind.right & Ind.bottom;
            y2Max          = this.optsNum.PhysArea.y2Max;    
            IntPathUpLow   = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
            
            IP                      = this.IC.SubShapePtsCart(this.RightCapillary.GetCartPts);           
            IntPath_Half_UpLow      = (this.RightCapillary.borderTop.IntSc - this.RightCapillary.borderBottom.IntSc)*IP;
            IntPath_Half_RightLeft  = (this.RightCapillary.borderRight.IntSc - this.RightCapillary.borderLeft.IntSc)*IP;
            
                        
            Diff           = this.IC.Diff;            
            phi_m          = this.optsPhys.phi_m; 
            M              = this.IC.M;
            F              = false(M,1);   
            T              = true(M,1);
            uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];
            
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            
            
            % Continuity   %[uv;phi;G]
            A_cont         = [diag(phi+phi_m)*Diff.div + [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)]...
                              diag(Diff.div*uv)+uvTdiag*Diff.grad,...
                              zeros(M)];
            v_cont         = Diff.div*(uv.*repmat(phi+phi_m,2,1));
                        
            
            % ************[uv;phi;G]
            ys             = fWP - Cn*Diff.Lap*phi;    
            ysP            = diag(fWPP) - Cn*Diff.Lap;
            Cmu            = -Diff.Lap*diag(phi+phi_m);
            
            Cuv            = Cak*(zeta + 4/3)*Diff.LapDiv; 
            
            A_cont(Ind.top|Ind.bottom,[T;T;F;F])  = Cuv(Ind.top|Ind.bottom,:);
            A_cont_phi                            = -Diff.Lap*diag(G) ...
                                                    + Diff.div*(diag(repmat(ys,2,1))*Diff.grad)...
                                                    + Diff.div*(diag(Diff.grad*phi)*repmat(ysP,2,1));
                                                    
            A_cont(Ind.top|Ind.bottom,[F;F;T;F])  = A_cont_phi(Ind.top|Ind.bottom,:);
            A_cont(Ind.top|Ind.bottom,[F;F;F;T])  = Cmu(Ind.top|Ind.bottom,:);
            
            v_cont(Ind.top|Ind.bottom)    =   Cmu(Ind.top|Ind.bottom,:)*G ...
                                            + Cuv(Ind.top|Ind.bottom,:)*uv ...
                                            + Diff.div(Ind.top|Ind.bottom,:)*(repmat(ys,2,1).*(Diff.grad*phi));
                                        
            % (BC1) y1 = +/- infinity
            A_cont(Ind.left|Ind.right,:)         = 0;
            A_cont(Ind.left|Ind.right,[F;F;F;T]) = Diff.Dy2(Ind.left|Ind.right,:);
            v_cont(Ind.left|Ind.right,:)         = Diff.Dy2(Ind.left|Ind.right,:)*G;                                        
            
            % ************            
            
            A_m11 = zeros(M,4*M);
            A_m11(:,[T;T;F;F]) = Cak*((zeta -2/3)*Diff.div + 2*[Diff.Dy1,zeros(M,M)]);
            A_m11(:,[F;F;T;F]) = diag(fWP - G) + Cn*( -diag(Diff.Dy1*phi)*Diff.Dy1 + diag(Diff.Dy2*phi).*Diff.Dy2);
            A_m11(:,[F;F;F;T]) = - diag(phi+phi_m);
            
            m11 =  Cak*((zeta -2/3)*Diff.div*uv + 2*Diff.Dy1*uv(1:end/2))...
                 + (fW - G.*(phi+phi_m) + Cn/2*( -(Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2 ));
            
            
            A_m12 = zeros(M,4*M);
            A_m12(:,[T;T;F;F]) = Cak*[Diff.Dy2 , Diff.Dy1];
            A_m12(:,[F;F;T;F]) = -Cn*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1);
            
            m12 = Cak*[Diff.Dy2 , Diff.Dy1]*uv - Cn*((Diff.Dy1*phi).*(Diff.Dy2*phi));
            
            
            A_cont(lbC,:) = IntPath_Half_UpLow*A_m12 +  IntPath_Half_RightLeft*A_m11;
            v_cont(lbC)   = IntPath_Half_UpLow*m12   +  IntPath_Half_RightLeft*m11;            
            
            % (BC2) [uv;phi;G]                   
            A_cont(rbC,:)   = IntPathUpLow*A_m12 + y2Max*(A_m11(rbC,:) - A_m11(lbC,:));
            v_cont(rbC)     = IntPathUpLow*m12 +y2Max*(m11(rbC) - m11(lbC));
                 
         end        
         function [v_mom,A_mom]   = Momentum(this,uv,phi,G)                
            %[uv;phi;G;p] 
            %
            %        A_mom          = [tauM,...
            %                          -[diag(Diff.Dy1*G);diag(Diff.Dy2*G)],...
            %                          - diag(repmat(phi+phi_m,2,1))*Diff.grad];
            % 
            %        v_mom         = tauM*uv - repmat(phi+phi_m,2,1).*(Diff.grad*G);

            Ind            = this.IC.Ind;
            Diff           = this.IC.Diff;            
            phi_m          = this.optsPhys.phi_m;             
            Cak            = this.optsPhys.Cak;    
            Cn             = this.optsPhys.Cn;
            zeta           = this.optsPhys.zeta;
            
            [fWP,fW,fWPP]  = DoublewellPotential(phi,Cn);
            tauM           = Cak*(Diff.LapVec + (zeta + 1/3)*Diff.gradDiv);            
            phiM2          = repmat(phi+phi_m,2,1);
           
            A_mom_phi      = -[diag(Diff.Dy1*G);diag(Diff.Dy2*G)]...
                             +diag(Diff.grad*phi)*repmat(diag(fWPP)-Cn*Diff.Lap,2,1)...
                             +diag(repmat(fWP - Cn*Diff.Lap*phi - G,2,1))*Diff.grad;       
            A_mom          = [tauM,...
                              A_mom_phi,...
                              - diag(phiM2)*Diff.grad - [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)]];
            v_mom         = tauM*uv - phiM2.*(Diff.grad*G) + ...
                             repmat(fWP - Cn*Diff.Lap*phi - G,2,1).*(Diff.grad*phi);
         end        
        
        function [v,A] = FullEqsBcs(this,z)
        end
    end
end