classdef DiffuseInterfaceSingleFluid < DiffuseInterface
    
    methods (Access = public)          
        SolveMovingContactLine(this,maxIterations)       
        [rho,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,rho,findTheta)
        [mu,uv,A,b,a]       = GetVelocityAndChemPot(this,rho,theta)
        [A,b]               = ContMom_DiffuseInterfaceSingleFluid(this,rho)                
        
        function p = GetPressure_from_ChemPotential(this,mu,rho_ig)
            Cn     = this.optsPhys.Cn;
            rho_m  = this.optsPhys.rho_m;
            %for bulk
            % 0 = W' - m
            % p = - W + mu*(rho + rho_m)
            
            fsolveOpts = optimset('Display','off');
            [rho,~,exitflag]   = fsolve(@f,rho_ig,fsolveOpts);
            if(exitflag < 1)
                cprintf('*r','No solution for density for given chemical potential found');
            end
            
            [~,W] = DoublewellPotential(rho,Cn);
            p     = - W + mu*(rho+rho_m);
            
            function y = f(rho)
                y = DoublewellPotential(rho,Cn) - mu;
            end
            
        end  
    end
end