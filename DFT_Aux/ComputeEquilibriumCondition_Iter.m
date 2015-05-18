function [sol] = ComputeEquilibriumCondition_Iter(params,misc)

    x_ig       = misc.x_ig;    
    Vext       = misc.Vext;
    VAdd       = misc.VAdd;
    Conv       = misc.Conv;
    IntMatrFex = misc.IntMatrFex;    
    nSpecies   = params.optsPhys.nSpecies;
    N          = size(Vext,1);
    
    if(~isfield(misc,'mark'))
        mark = true(size(Vext,1),1);
    else
        mark = misc.mark;
    end
    
    kBT    = params.optsPhys.kBT;    
    
	getFex    = str2func(['Fex_',params.FexNum.Fex]);    
    if(isfield(params.optsPhys,'sigmaS'))        
        R      = params.optsPhys.sigmaS/2;
    else
        R = [];
    end
    
    fsolveOpts       = optimset('TolFun',1e-8,'TolX',1e-8);
    
    if(isfield(params.optsPhys,'mu_sat') && isfield(params.optsPhys,'Dmu'))
       params.optsPhys.mu =  params.optsPhys.mu_sat + params.optsPhys.Dmu;       
    end
    
    if(isfield(params.optsPhys,'mu'))   
        mu                      = params.optsPhys.mu;
        [x_ic,~,flag_fsolve]    = fsolve(@fs,x_ig(mark),fsolveOpts);
    else        
        Int                     = misc.Int;
        nParticlesS             = params.optsPhys.nParticlesS;
        [x_ic,errorHistory1]     = NewtonMethod(x_ig([true;mark],:),@fs_canonical,1,100,0.05);
        [x_ic,errorHistory2]     = NewtonMethod(x_ic,@fs_canonical,1e-10,100,0.5);
        %[x_ic,~,flag_fsolve]    = fsolve(@fs_canonical,x_ig([true;mark],:),fsolveOpts);
        
        if(isempty(x_ic))
            disp('error');
        else
            mu    = x_ic(1,:);
            x_ic  = x_ic(2:end,:);        
        end
    end
    
    x_ic_full(mark,:)  = x_ic;
    x_ic_full(~mark,:) = x_ig(~mark,:);
    
    rho = exp((x_ic_full-Vext)/kBT);   
    
    sol   = v2struct(rho,mu);    
    sol.x = x_ic_full;
    
    function mu_sRel = fs(xm)                
        mu_sRel = GetExcessChemPotentialPart(xm,mu);%./exp((xm-Vext(mark,:))/kBT);
    end

    function [y,J] = fs_canonical(x)
        mu_s         = x(1,:);
        x            = x(2:end,:);       
        
        [y,J]        = GetExcessChemPotentialPart(x,mu_s);%./exp((x-Vext(mark,:))/kBT);
        
        xf(mark,:)   = x;
        xf(~mark,:)  = x_ig(~mark,:);        
        rho_full     = exp((xf-Vext)/kBT);
        y            = [Int*rho_full - nParticlesS';y];
        y            = y(:);   
        
        J            = [0,Int*diag(rho_full)/kBT;J];
    end

    function [mu_s,J_s] = GetExcessChemPotentialPart(xm,mu)        
        x(mark,:)  = xm;
        x(~mark,:) = x_ig(~mark,:);        
        [mu_s,J_s] = GetExcessChemPotential(x,mu);
        mu_s       = mu_s(mark,:);
        J_s        = J_s(mark,[true;mark]);
    end   
    function [mu_s,J_s] = GetExcessChemPotential(x,mu)
        rho_s            = exp((x-Vext)/kBT);
        [mu_HS,~,~,J_HS] = getFex(rho_s,IntMatrFex,kBT,R);
        if(size(J_HS,2) == nSpecies)
            J_HS = diag(J_HS);
        end        
        [mu_attr] = Fex_Meanfield(rho_s,Conv,kBT);
                            
        %mu_s = mu_HS + mu_attr + x + VAdd;
        %J_s  = (J_HS + Conv.Conv)*diag(rho_s)/kBT + eye(N);
        
        mu_s = mu_HS + Conv.Conv*rho_s + x + VAdd;
        J_s  = (J_HS + Conv.Conv)*diag(rho_s)/kBT + eye(N);
        
        mu_s = mu_s - mu;
        J_s  = [-ones(N,1),J_s];
        %for iSpecies=1:nSpecies
        %   mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu(iSpecies);
        %end
        
    end   

end