classdef DiffuseInterfaceBinaryFluid < DiffuseInterface
    
   properties (Access = public)
       p %pressure       
   end
       
	methods (Access = public) 
       function this = DiffuseInterfaceBinaryFluid(config)           
             this@DiffuseInterface(config);
       end
       
       [rho,theta] = GetEquilibriumDensity(this,mu,theta,rho,findTheta)
       [A,b] = ContMom_DiffuseInterfaceBinaryFluisFluid(this,phi_s)
   end
end