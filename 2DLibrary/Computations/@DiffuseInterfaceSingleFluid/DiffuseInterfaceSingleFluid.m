classdef DiffuseInterfaceSingleFluid < DiffuseInterface
    
    methods (Access = public) 
         function this = DiffuseInterfaceSingleFluid(config)           
             this@DiffuseInterface(config);
         end
        
        SolveMovingContactLine(this,maxIterations)       
        [phi,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,phi,findTheta)
        [mu,uv,A,b,a]       = GetVelocityAndChemPot(this,phi,theta)
        [A,b]               = [A,b] = ContinuumMomentumEqs(this,phi)
                        
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
        
        function [bulkError,bulkAverageError] = DisplayFullError(this,phi,uv)            
            M        = this.IC.M;
            
            [Af,bf]  = ContinuumMomentumEqs(this,phi);                         
            mu_s     = GetMu(this,phi);
            
            error    = Af*[mu_s;uv] - bf;
            PrintErrorPos(error(1:M),'continuity equation',this.IC.Pts);            
            PrintErrorPos(error(1+M:2*M),'y1-momentum equation',this.IC.Pts);            
            PrintErrorPos(error(2*M+1:end),'y2-momentum equation',this.IC.Pts);                        
            
            PrintErrorPos(error(repmat(~this.IC.Ind.bound,2,1)),'bulk continuity and momentum equations');               
            
            bulkError        = max(abs(error(repmat(~this.IC.Ind.bound,2,1))));
            bulkAverageError = mean(abs(error(repmat(~this.IC.Ind.bound,2,1))));
        end
    end
end