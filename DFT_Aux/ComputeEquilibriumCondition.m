function [x_ic] = ComputeEquilibriumCondition(params,misc)


    x_ig = misc.x_ig;
    mark = misc.mark;
    Vext = misc.Vext;
    VAdd = misc.VAdd;
    Conv = misc.Conv;
    IntMatrFex = misc.IntMatrFex;
    
    kBT    = params.optsPhys.kBT;
    mu_sat = params.optsPhys.mu_sat;    
    Dmu    = params.optsPhys.Dmu;    
    R      = params.optsPhys.sigmaS/2;     
    
	getFex    = str2func(['Fex_',misc.FexNum.Fex]);
    
    fsolveOpts       = optimset('TolFun',1e-10,'TolX',1e-10);
    
    %[x_ic,~,flag_fsolve]    = fsolve(@GetExcessChemPotential,x_ig(mark),fsolveOpts);
    [x_ic,~,flag_fsolve]    = fsolve(@fs,x_ig(mark),fsolveOpts);
    if(flag_fsolve~=1)
        x_ic = 0;
    end
    
    function mu_sRel = fs(xm)
        mu_sRel = GetExcessChemPotential(xm)./exp((xm-Vext(mark))/kBT);
    end
    function mu_s = GetExcessChemPotential(xm)        
        x(mark)  = xm;
        x(~mark) = x_ig(~mark);
        x        = x';
        mu_s = GetExcessChemPotentialFull(x);
        mu_s = mu_s(mark);
    end   
    function mu_s = GetExcessChemPotentialFull(x)                
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R);
                       
%        for iSpecies=1:1%nSpecies
%           mu_s(:,iSpecies) = mu_s(:,iSpecies) - (mu_sat + Dmu);
        %end
        mu_s = mu_s - (mu_sat + Dmu);
        mu_s = mu_s + x + Conv*rho_s + VAdd;           
    end   
end