function [sol] = ComputeEquilibriumCondition(params,misc)

    x_ig       = misc.x_ig;    
    Vext       = misc.Vext;
    VAdd       = misc.VAdd;
    Conv       = misc.Conv;
    Int        = misc.Int;
    IntMatrFex = misc.IntMatrFex;    
    nSpecies   = params.optsPhys.nSpecies;
    N          = size(Vext,1);            
    
    
    if(~isfield(params,'solver'))
        params.solver = 'fsolve';
    end
    
    if(~isfield(misc,'mark'))
        mark = true(size(Vext,1),1);
    else
        mark = misc.mark;
    end    
    
    kBT    = params.optsPhys.kBT;    
    
    if(isfield(params,'optsNum'))
        FexNum = params.optsNum.FexNum;
    else
        FexNum = params.FexNum;
    end
    
	getFex    = str2func(['Fex_',FexNum.Fex]);    
    if(isfield(params.optsPhys,'sigmaS'))        
        R      = params.optsPhys.sigmaS/2;
    else
        R = [];
    end        
    
    if(isfield(params.optsPhys,'mu_sat') && isfield(params.optsPhys,'Dmu'))
       params.optsPhys.mu =  params.optsPhys.mu_sat + params.optsPhys.Dmu;       
    end
    
    if(strcmp(FexNum.Fex,'FMTRosenfeld_3DFluid') && strcmp(params.solver,'Newton'))
        fullInput = true; %case only implemented for 1 species
    else
        fullInput = false;
    end

    if(fullInput)
        if(isfield(params.optsPhys,'mu'))
            rho_ig   = exp((x_ig-Vext)/kBT);
        else
            rho_ig   = exp((x_ig(2:end)-Vext)/kBT);
        end
        x_ig     = [x_ig;     
                    IntMatrFex(1).AD.n2 * rho_ig;...
                    IntMatrFex(1).AD.n3 * rho_ig;...
                    IntMatrFex(1).AD.n2_v_1 * rho_ig;...
                    IntMatrFex(1).AD.n2_v_2 * rho_ig];       
        
        N_AD     = size(IntMatrFex(1).AD.n2,1); %size(n2_i,1);
        NFull    = N+4*N_AD;         
        if(isfield(params,'maxComp_y2'))
            markFull = [(misc.PtsCart.y2_kv <= params.maxComp_y2);...
                         repmat(misc.AD.PtsCart.y2_kv <= params.maxComp_y2+R,4,1)];   
        else
            markFull = true(N+4*N_AD,1);
        end
    else
        markFull = mark;
        NFull    = N;
    end            
    
    if(isfield(params.optsPhys,'mu'))
        %x_ig_n      = zeros(NFull,nSpecies);
        %x_ig_n(1:N) = x_ig;%(2:end);
        x_ig_n                 = x_ig;%(2:end);
        mu                     = params.optsPhys.mu;
        
        %Picard-Iteration        
        if(strcmp(params.solver,'Picard'))        
            markFull = mark;
            it   = 1; err  = 1;
            x_ig_n = x_ig_n(1:N);
            %x_ic = x_ic(markFull);            
            x_ic = x_ig_n(markFull);
            
            no = 0;
            while((err > 0.1) && (it < 1000))
                [dx,err]   = GetdX_nP1(x_ic);
                x_ic       = x_ic + 0.01*dx;
                
                for ih = 1:no
                        fprintf('\b');
                end 
                no = fprintf(['Iteration ',num2str(it),': ',num2str(err)]);
                it = it+1;
            end
            fprintf('\n'); no = 0;

            while((err > 1e-4) && (it < 1000))
                [dx,err]   = GetdX_nP1(x_ic);
                x_ic       = x_ic + 0.1*dx;
                for ih = 1:no
                        fprintf('\b');
                end 
                no = fprintf(['Iteration ',num2str(it),': ',num2str(err)]);
                it = it+1;
            end
            fprintf('\n'); no = 0;
             while((err > 1e-10) && (it < 50000))
                [dx,err]   = GetdX_nP1(x_ic);
                x_ic       = x_ic + 0.2*dx;
                for ih = 1:no
                        fprintf('\b');
                end 
                no = fprintf(['Iteration ',num2str(it),': ',num2str(err)]);
                it = it+1;
             end                  
             fprintf('\n'); no = 0;
        elseif(strcmp(params.solver,'Newton'))
            [x_ic,errorHistory1]    = NewtonMethod(x_ig_n(markFull),@fs,1e-10,20,0.7,{'returnLastIteration'});
            [x_ic,errorHistory2]    = NewtonMethod(x_ic,@fs,1e-10,10000,1,{'returnLastIteration'});
            errorHistory            = [errorHistory1 errorHistory2];
        else
            fsolveOpts             = optimset('TolFun',1e-8,'TolX',1e-8);
            [x_ic,~,flag_fsolve]   = fsolve(@fs,x_ig(mark),fsolveOpts);
        end        
        
    else                
        nParticlesS             = params.optsPhys.nParticlesS;     
                
        if(strcmp(params.solver,'Newton'))    
            x_ig_n = zeros(NFull+1,nSpecies);
            %[x_ic,errorHistory1]    = NewtonMethod(x_ig_n(:),@fs_canonical,1,100,0.7);
            [x_ic,errorHistory1]    = NewtonMethod(x_ig_n(:),@fs_canonical,1,100,0.2);
            %[x_ic,errorHistory2]    = NewtonMethod(x_ic,@fs_canonical,1e-10,200,1,{'returnLastIteration'});
            [x_ic,errorHistory2]    = NewtonMethod(x_ic,@fs_canonical,1e-10,200,1,{'returnLastIteration'});
            errorHistory            = [errorHistory1 errorHistory2];
        else
            x_ig_n                  = x_ig([true;mark],:);
            fsolveOpts              = optimset('TolFun',1e-8,'TolX',1e-8);
   
            [x_ic,~,flag_fsolve]    = fsolve(@fs_canonical,x_ig_n(:),fsolveOpts);
            %[x_ic,~,flag_fsolve]    = fsolve(@fs_canonical,x_ig([true;mark],:),fsolveOpts);
        end        
        
        if(isempty(x_ic))
            disp('error');
        else
            mu    = (x_ic(1:nSpecies))';
            x_ic  = x_ic(nSpecies+1:end);            
        end
    end
    
    x_ic      = reshape(x_ic,[],nSpecies);        
    x_ic      = GetFullX(x_ic);%    
    rho       = exp((x_ic(1:N,:)-Vext)/kBT);
    
    %x_ic_full(mark,:)  = x_ic;
    %x_ic_full(~mark,:) = x_ig(~mark,:);
%    
%    rho = exp((x_ic_full-Vext)/kBT);   
%    
    sol   = v2struct(rho,mu);    
    sol.x = x_ic(1:N,:);
    
    function [mu_sRel,J] = fs(xm)      
        if(nargout >= 2)
            [mu_sRel,J] = GetExcessChemPotentialPart(xm,mu);%./exp((xm-Vext(mark,:))/kBT);
            J           = J(:,2:end);
        else
            [mu_sRel] = GetExcessChemPotentialPart(xm,mu);%./exp((xm-Vext(mark,:))/kBT);
            J         = [];
        end        
    end
    function [y,J] = fs_canonical(x)
        mu_s         = (x(1:nSpecies))';
        x            = x(nSpecies+1:end);       
        x            = reshape(x,NFull,nSpecies);
        
        if(nargout >= 2)
            [y,J]    = GetExcessChemPotentialPart(x,mu_s);%./exp((x-Vext(mark,:))/kBT);
        else
            y        = GetExcessChemPotentialPart(x,mu_s);%./exp((x-Vext(mark,:))/kBT);
        end
        
        xf           = GetFullX(x);
        rho_full     = exp((xf(1:N,:)-Vext)/kBT);
        
        %Add mass constraint
        yP           = Int*rho_full - nParticlesS';        
        y            = [yP(:);y(:)];        
        
        if(nargout >= 2)
            JP = zeros(nSpecies,nSpecies*NFull);
            for i = 1:nSpecies
                JP(i,1+(i-1)*N:i*N) = Int*diag(rho_full(:,i))/kBT;
            end        
            J            = [zeros(nSpecies),JP;J];
        end
    end
    function [mu_s,J_s] = GetExcessChemPotentialPart(xm,mu)        
        x          = GetFullX(xm);
        if(nargout >= 2)
            [mu_s,J_s] = GetExcessChemPotential(x,mu);
        else
            [mu_s] = GetExcessChemPotential(x,mu);
            J_s    = [];
        end
        
        markF      = repmat(markFull,nSpecies,1);
        mu_s       = mu_s(markFull,:);
        %mu_s       = mu_s(:);
        if(~isempty(J_s))
            J_s        = J_s(markF,[true(nSpecies,1);markF]);
        end
    end   
    function [mu_s,J_s] = GetExcessChemPotential(x,mu)
        xR              = (x(1:N,:));        
        rho_s           = exp((xR-Vext)/kBT);
        x(1:N,:)        = rho_s;                
        %rho_s            = exp((x-Vext)/kBT);
        if(nargout >= 2)
            [mu_HS,~,~,J_HS] = getFex(x,IntMatrFex,kBT,R);
            if(size(J_HS,2) == nSpecies)
                J_HS = diag(J_HS(:));
            end
        else
            [mu_HS]          = getFex(x,IntMatrFex,kBT,R);
            J_HS = [];
        end

        [mu_attr,J_attr] = Fex_Meanfield(rho_s,Conv,kBT);
                            
        mu_s             = mu_HS;
        mu_s(1:N,:)      = mu_HS(1:N,:) + mu_attr + xR + VAdd;        
        
        for iSpecies=1:nSpecies
           mu_s(1:N,iSpecies) = mu_s(1:N,iSpecies) - mu(iSpecies);                      
        end        
        
        J_s                     = J_HS;
        if(~isempty(J_s))
            markRhos = (1:N*nSpecies);%repmat(markRhos,nSpecies,1);

            J_s                     = J_HS;
            J_s(markRhos,markRhos)  = (J_HS(markRhos,markRhos) + J_attr)*diag(rho_s(:))/kBT + eye(N*nSpecies);                   

            J_s = [zeros(size(J_s,1),nSpecies),J_s];
            for iSpecies=1:nSpecies              
               J_s((1+(iSpecies-1)*N):(iSpecies*N),iSpecies) = -ones(N,1);
            end
        end
    end   
    function [dx,err] = GetdX_nP1(x)
        x       = GetFullX(x);
        dx      = -GetExcessChemPotential(x,mu);    
        dx      = dx(markFull);
        err     = max(abs(dx));
    end

    function xFull = GetFullX(x)
        xFull(markFull,:)  = x;
        xFull(~markFull,:) = x_ig(~markFull,:);        
    end

end