function [sol] = ComputeEquilibriumCondition_Iter(params,misc)

    x_ig       = misc.x_ig;    
    Vext       = misc.Vext;
    VAdd       = misc.VAdd;
    Conv       = misc.Conv;
    Int        = misc.Int;
    IntMatrFex = misc.IntMatrFex;    
    nSpecies   = params.optsPhys.nSpecies;
    N          = size(Vext,1);
    
    if(~isfield(misc,'mark'))
        mark = true(size(Vext,1),1);
    else
        mark = misc.mark;
    end
    markF = repmat(mark,nSpecies,1);
    
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
        [x_ic,errorHistory]    = NewtonMethod(x_ig(mark),@fs,1e-10,100,1);
        %[x_ic,errorHistory2]     = NewtonMethod(x_ic,@fs,1e-10,100,1);
        %[x_ic,~,flag_fsolve]    = fsolve(@fs,x_ig(mark),fsolveOpts);
    else                
        nParticlesS             = params.optsPhys.nParticlesS;
        x_ig_n                  = x_ig([true;mark],:);
        [x_ic,errorHistory1]    = NewtonMethod(x_ig_n(:),@fs_canonical,10,100,0.7);
        [x_ic,errorHistory2]     = NewtonMethod(x_ic,@fs_canonical,1e-10,100,1);
        %[x_ic,~,flag_fsolve]    = fsolve(@fs_canonical,x_ig([true;mark],:),fsolveOpts);
        errorHistory = [errorHistory1 errorHistory2];
        
        if(isempty(x_ic))
            disp('error');
        else
            mu    = (x_ic(1:nSpecies))';
            x_ic  = x_ic(nSpecies+1:end);
            x_ic  = reshape(x_ic,N,nSpecies);
        end
    end
    
    x_ic_full(mark,:)  = x_ic;
    x_ic_full(~mark,:) = x_ig(~mark,:);
    
    rho = exp((x_ic_full-Vext)/kBT);   
    
    sol   = v2struct(rho,mu);    
    sol.x = x_ic_full;
    
    function [mu_sRel,J] = fs(xm)                
        [mu_sRel,J] = GetExcessChemPotentialPart(xm,mu);%./exp((xm-Vext(mark,:))/kBT);
        J           = J(:,2:end);
    end

    function [y,J] = fs_canonical(x)
        mu_s         = (x(1:nSpecies))';
        x            = x(nSpecies+1:end);       
        x            = reshape(x,N,nSpecies);
        
        [y,J]        = GetExcessChemPotentialPart(x,mu_s);%./exp((x-Vext(mark,:))/kBT);
        
        xf(mark,:)   = x;
        xf(~mark,:)  = x_ig(~mark,:);        
        rho_full     = exp((xf-Vext)/kBT);
        
        %Add mass constraint
        yP           = Int*rho_full - nParticlesS';        
        y            = [yP(:);y(:)];        
        
        JP = zeros(nSpecies,nSpecies*N);
        for i = 1:nSpecies
            JP(i,1+(i-1)*N:i*N) = Int*diag(rho_full(:,i))/kBT;
        end        
        J            = [zeros(nSpecies),JP;J];
    end

    function [mu_s,J_s] = GetExcessChemPotentialPart(xm,mu)        
        x(mark,:)  = xm;
        x(~mark,:) = x_ig(~mark,:);        
        [mu_s,J_s] = GetExcessChemPotential(x,mu);
        mu_s       = mu_s(mark,:);
        J_s        = J_s(markF,[true(nSpecies,1);markF]);
    end   
    function [mu_s,J_s] = GetExcessChemPotential(x,mu)
        rho_s            = exp((x-Vext)/kBT);
        [mu_HS,~,~,J_HS] = getFex(rho_s,IntMatrFex,kBT,R);
        if(size(J_HS,2) == nSpecies)
            J_HS = diag(J_HS(:));
        end        
        [mu_attr,J_attr] = Fex_Meanfield(rho_s,Conv,kBT);
                            
        mu_s = mu_HS + mu_attr + x + VAdd;
        J_s  = (J_HS + J_attr)*diag(rho_s(:))/kBT + eye(N*nSpecies);
        
                        
        J_s = [zeros(nSpecies*N,nSpecies),J_s];
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies)                     = mu_s(:,iSpecies) - mu(iSpecies);           
           J_s((1+(iSpecies-1)*N):(iSpecies*N),iSpecies) = -ones(N,1);
        end
    end   

end