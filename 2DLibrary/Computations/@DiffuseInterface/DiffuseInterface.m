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
            this.IC.SetUpBorders(300);            
        end
        
%         function [z,dz,ddz] = W(this,rho)        
%             Cn  = this.optsPhys.Cn;
%             z   = (1-rho.^2).^2/(2*Cn);
%             dz  = -2*rho.*(1-rho.^2)/Cn;
%             ddz = 2*(3*rho.^2-1)/Cn;
%         end
        
        function [A,b] = CahnHilliard_StressTensorIJ(this,rho,i,j)
            % get matrices for
            % T = Ca*( eta*(grad(u) + grad(u)^T) + (zeta - 2/3*eta) div(u)*I ) +...
            %           + (W(rho) + Cn/2*|grad(rho)|^2 - mu*(rho+rho_m))*I - Cn*(grad(rho) X grad(rho))
            %   = A(rho)*[mu;uv] + b(rho)
            
            Cn    = this.optsPhys.Cn;
            Cak   = this.optsPhys.Cak;
            eta   = this.optsPhys.eta;            
            zeta  = this.optsPhys.zeta;
            rho_m = this.optsPhys.rho_m;
            Diff  = this.IC.Diff;
            M     = this.IC.M;
            
            [~,W]   = DoublewellPotential(rho,Cn);
            bDiag   = W + Cn/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
            %bDiag   = W(rho) + 1/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
            if(i == 1 && j == 1)
                Amu     = -diag(rho+rho_m);
                Auv     = Cak*(eta*[2*Diff.Dy1 , zeros(M)] + (zeta - 2/3*eta)*Diff.div);
                b       = bDiag - Cn*(Diff.Dy1*rho).*(Diff.Dy1*rho);
            elseif((i==2 && j == 1)  || (i==1 && j == 2))
                Amu     = zeros(M);
                Auv     = Cak*[Diff.Dy2 , Diff.Dy1]*eta;
                b       = - Cn*(Diff.Dy1*rho).*(Diff.Dy2*rho);
            elseif(i==2 && j == 2)
                Amu     = -diag(rho+rho_m);
                Auv     = Cak*(eta*[zeros(M) , 2*Diff.Dy2] + (zeta - 2/3)*Diff.div);
                b       = bDiag - Cn*(Diff.Dy2*rho).*(Diff.Dy2*rho);
            end

            A = [Amu Auv];

        end                     
        
        function [A,b] = CahnHilliard_DivergenceOfStressTensor(this,rho)

            Diff      = this.IC.Diff;    
            
            [A11,b11] = CahnHilliard_StressTensorIJ(this,rho,1,1); 
            [A12,b12] = CahnHilliard_StressTensorIJ(this,rho,1,2); 
            [A21,b21] = CahnHilliard_StressTensorIJ(this,rho,2,1); 
            [A22,b22] = CahnHilliard_StressTensorIJ(this,rho,2,2); 

            A1        = Diff.Dy1 * A11  + Diff.Dy2*A21;
            A2        = Diff.Dy1 * A12  + Diff.Dy2*A22;

            b1        = Diff.Dy1 * b11  + Diff.Dy2*b21;
            b2        = Diff.Dy1 * b12  + Diff.Dy2*b22;

            A   = [A1;A2];
            b   = [b1;b2];      

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
        function [A,b] = ContMom_DiffuseInterfaceSingleFluid(this,rho)
        %Continuitiy: div(rho*uv) = 0
        %Momentum: - (rho+rho_m)*grad(mu) +
        %          ... +grad(rho)*(W'(rho) - mu - Cn*Lap(rho) )
        %          + Cak*(eta*Lap(uv) + (zeta + eta/3)*grad(div(uv)) ) +...

        %
        % A*[mu;uv] = b corresponds to momentum and continuity Eq. for given rho               
            rho_m          = this.optsPhys.rho_m;
            Cn             = this.optsPhys.Cn;
            Cak            = this.optsPhys.Cak;
            eta            = this.optsPhys.eta;
            zeta           = this.optsPhys.zeta;
            Diff           = this.IC.Diff;
            M              = this.IC.M;
            
            rho_f          = rho + rho_m;
            rho_f2         = repmat(rho_f,2,1);
            gradRho_T      = Diff.grad*rho_f;

            aT             = zeros(M,2*M);
            aT(:,1:M)      = diag(gradRho_T(1:M));
            aT(:,1+M:end)  = diag(gradRho_T(1+M:end));

            A_cont_mu      = zeros(M);
            A_cont_uv      = aT + diag(rho_f)*Diff.div; 

            A_mom_mu       = -diag(rho_f2)*Diff.grad - [diag(Diff.Dy1*rho);diag(Diff.Dy2*rho)];
            A_mom_uv       = Cak*(eta*Diff.LapVec + (zeta + eta/3)*Diff.gradDiv);

            A_cont         = [A_cont_mu,A_cont_uv];
            A_mom          = [A_mom_mu, A_mom_uv];
            A              = [A_cont;A_mom];  

            b              = zeros(3*M,1);            
            ys             = DoublewellPotential(rho,Cn) - Cn*(Diff.Lap*rho);

            b(1+M:end)     = - repmat(ys,2,1).*(Diff.grad*rho); 

        end        
        %**************************************
        %Plot Functions
        %**************************************
        function PlotSeppecherSolution(this,D_B,theta,rho)
            
            UWall   = this.optsPhys.UWall;
            D_A     = this.optsPhys.D_A;
            PtsCart = this.IC.GetCartPts();
            y2Max   = this.optsNum.PhysArea.y2Max;
            Ind     = this.IC.Ind;
                        
            u_flow = GetSeppecherSolutionCart(PtsCart,UWall,D_A,D_B,theta);                        
            
            figure('Position',[0 0 1800 600],'Color','white');
            subplot(1,2,1);  this.IC.doPlotsStreamlines(u_flow,Ind.top); %IC.doPlotsFlux(u_flow)                        
            
            if(nargin >= 4)
                subplot(1,2,2);  this.IC.doPlotFLine([-100,100],[y2Max,y2Max],rho.*u_flow(end/2+1:end),'CART'); 
            end
            title('Check accuracy of map');
        end
        
        function PlotMu_and_U(this,mu,uv)
            
            Ind = this.IC.Ind;
            
            figure('Position',[0 0 1800 600],'color','white');
            subplot(1,2,1); this.IC.doPlots(mu);
            subplot(1,2,2); this.IC.doPlotsStreamlines(uv,Ind.top|Ind.right|Ind.left); %IC.doPlotsFlux(u_flow)(mu);
        end
           
        %**********************************
        %External files:
        %**********************************
        rho = GetEquilibriumDensity(this,theta,nParticles,rho)
        D_B = SetD_B(this,theta,rho,initialGuessDB)
        [mu,uv,A,b] = GetVelocityAndChemPot(this,rho,D_B,theta)
   end
    
end