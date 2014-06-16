classdef DiffuseInterface < handle
    
    
   properties (Access = public)
              
       IC       
       
       optsNum
       optsPhys
                                    
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
        
        function [z,dz,ddz] = W(this,rho)        
            Cn  = this.Cn;
            z   = (1-rho.^2).^2/(2*Cn);
            dz  = -2*rho.*(1-rho.^2)/Cn;
            ddz = 2*(3*rho.^2-1)/Cn;
        end
        
        function [A,b] = CahnHilliard_StressTensorIJ(rho,i,j)
            % get matrices for
            % T = Ca*( eta*(grad(u) + grad(u)^T) + (zeta - 2/3*eta) div(u)*I ) +...
            %           + (W(rho) + Cn/2*|grad(rho)|^2 - mu*(rho+rho_m))*I - Cn*(grad(rho) X grad(rho))
            %   = A(rho)*[mu;uv] + b(rho)
            

            bDiag   = W(rho) + Cn/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
            %bDiag   = W(rho) + 1/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);
            if(i == 1 && j == 1)
                Amu     = -diag(rho+rho_m);
                Auv     = Ca*(eta*[2*Diff.Dy1 , zeros(M)] + (zeta - 2/3*eta)*Diff.div);
                b       = bDiag - Cn*(Diff.Dy1*rho).*(Diff.Dy1*rho);
            elseif((i==2 && j == 1)  || (i==1 && j == 2))
                Amu     = zeros(M);
                Auv     = Ca*[Diff.Dy2 , Diff.Dy1]*eta;
                b       = - Cn*(Diff.Dy1*rho).*(Diff.Dy2*rho);            
            elseif(i==2 && j == 2)
                Amu     = -diag(rho+rho_m);
                Auv     = Ca*(eta*[zeros(M) , 2*Diff.Dy2] + (zeta - 2/3)*Diff.div);
                b       = bDiag - Cn*(Diff.Dy2*rho).*(Diff.Dy2*rho);
            end

            A = [Amu Auv];

        end
            
   end
    
end