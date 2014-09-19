classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p %pressure       
   end
       
	methods (Access = public) 
       function this = DiffuseInterfaceBinaryFluid(config)           
             this@DiffuseInterface(config);
       end
       
       [rho,theta] = GetEquilibriumDensity(this,theta,phi,uv)
       [A,b] = ContinuumMomentumEqs(this,phi)
       SolveMovingContactLine(this,noIterations)  
       [p,uv,A,b,a] = GetVelocityAndChemPot(this,phi,theta)    
       
       function [bulkError,bulkAverageError] = DisplayFullError(this)
            M        = this.IC.M;
            m        = this.optsPhys.mobility;
            Cn       = this.optsPhys.Cn;
            Diff     = this.IC.Diff;
            
            uv       = this.uv;
            phi      = this.phi;
            p        = this.p;
            
            [Af,bf]  = ContinuumMomentumEqs(this,phi);                                     
            Convect    = [diag(uv(1:end/2)),diag(uv(1+end/2:end))]*this.IC.Diff.grad/m;
            
            error    = Af*[p;uv] - bf;
            
            absGradRho  = ( (Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2 );
            errorPhi    =  - Convect*phi ...
                            + 12*Cn*phi.*absGradRho ...
                            + 2*Cn*(3*phi.^2 -1).*(Diff.Lap*phi)...
                            - Cn*Diff.Lap2*phi;
            %Convect*phi - this.IC.Diff.Lap*GetMu(this);
            
            PrintErrorPos(error(1:M),'continuity equation',this.IC.Pts);            
            PrintErrorPos(error(1+M:2*M),'y1-momentum equation',this.IC.Pts);            
            PrintErrorPos(error(2*M+1:end),'y2-momentum equation',this.IC.Pts);                        
            PrintErrorPos(errorPhi,'Phasefield equation',this.IC.Pts);
            
            PrintErrorPos(error(repmat(~this.IC.Ind.bound,3,1)),'bulk continuity and momentum equations');               
            PrintErrorPos(errorPhi(~this.IC.Ind.bound),'bulk phasefield equations');               
            
            bulkError        = max(abs(error(repmat(~this.IC.Ind.bound,2,1))));
            bulkAverageError = mean(abs(error(repmat(~this.IC.Ind.bound,2,1))));
        end
   end
end