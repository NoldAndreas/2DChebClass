function ComputeDynamicsInertia(this,x_ic,mu)

%% Solves
% 
% $$\frac{\partial \varrho}{\partial t} = \nabla \cdot ( \varrho \nabla \mu )$$
% 
% in computational variable x: $\varrho = e^{(x-Vext)/(k_B T)}$
% 
% $$\frac{\partial x}{\partial t} = k_B T \Delta \mu + (\nabla x - \nabla V)\cdot \nabla \mu$$
% 

    % Initialization
    optsPhys    = this.optsPhys;
    optsNum     = this.optsNum;
    
    kBT         = optsPhys.kBT;
    if(isfield(optsPhys,'sigmaS'))    
        R       = optsPhys.sigmaS/2;
    else
        R       = [];
    end
    gammaS      = optsPhys.gammaS;
    mS          = optsPhys.mS;
    Diff        = this.IDC.Diff;
    plotTimes   = this.optsNum.plotTimes;
    nSpecies    = this.optsPhys.nSpecies;
    M           = this.IDC.M;
    Vext        = this.Vext;
    Vext_grad   = this.Vext_grad;
    IntMatrHI   = this.IntMatrHI;
    IntMatrFex  = this.IntMatrFex;
    Int_of_path = this.Int_of_path;
    Conv        = this.IntMatrV2;
    Ind         = this.IDC.Ind;        
    getFex      = str2func(['Fex_',optsNum.FexNum.Fex]);    
    doHI        = this.doHI;
    
    markVinf    = (Vext == inf);
    if(strcmp(this.IDC.polar,'polar'))
        polarShape = true;
    else
        polarShape = false;
    end
    
    subArea     = this.subArea;
    
	I           = eye(M);  
    eyes        = [I I];
    
    PtsCart = this.IDC.GetCartPts();
	y1S     = repmat(PtsCart.y1_kv,1,nSpecies); 
    y2S     = repmat(PtsCart.y2_kv,1,nSpecies);

    ythS    = repmat(this.IDC.Pts.y2_kv,1,nSpecies);
    
    if(nargin < 2)
        x_ic = this.x_eq;
        mu   = this.mu;
    end
        
    %tic   
    
    fprintf(1,'Computing dynamics ...'); 

    if(isfield(optsNum,'PlotArea'))
        optsNumT = rmfield(optsNum,'PlotArea');
    end
    [this.dynamicsResult,recEq,paramsEq] = DataStorage('Dynamics',...
                            @ComputeDDFTDynamics,v2struct(optsNumT,optsPhys),[]); %true      
                     
    function data = ComputeDDFTDynamics(params,misc)        
        
        finite1 = (Ind.normalFinite*[ones(M,1);zeros(M,1)]~=0);
        finite2 = (Ind.normalFinite*[zeros(M,1);ones(M,1)]~=0);

        mMx              = ones(M,1);        
        mMx(markVinf(:,1)) = 0;  % assumes infinite potential same for each species
        mMuv              = ones(2*M,1);
        %mMuv([Ind.finite;Ind.finite]) = 0;
        mMuv([finite1;Ind.finite]) = 0;
        
        mM = [mMx;mMuv];
        
        mM = repmat(mM,nSpecies,1);

        x_ic = [x_ic;zeros(2*M,nSpecies)]; % padd with zero velocity
        
        opts    = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([ones(nSpecies,1);mM]));    
        [~,X_t] = ode15s(@dx_dt,plotTimes,[zeros(nSpecies,1);x_ic(:)],opts);   
                       
        nPlots    = length(plotTimes);

        accFlux   = X_t(:,1:nSpecies);
        X_t       = X_t(:,nSpecies+1:end)';

        Y_t       = X_t(1:M,:);
        UV_t      = X_t(M+1:3*M,:);
        
        rho_t     = exp((Y_t-Vext(:)*ones(1,nPlots))/kBT);

        X_t       = reshape(Y_t,M,nSpecies,nPlots);
        UV_t      = reshape(UV_t,M,nSpecies,nPlots);
        rho_t     = reshape(rho_t,M,nSpecies,nPlots);
        flux_t    = zeros(2*M,nSpecies,nPlots);
        V_t       = zeros(M,nSpecies,nPlots);   

        for i = 1:length(plotTimes)
            flux_t(:,:,i) = GetFlux(X_t(:,:,i),plotTimes(i));           
            V_t(:,:,i)    = Vext + getVAdd(y1S,y2S,plotTimes(i),optsPhys.V1);
        end

        data       = v2struct(IntMatrFex,X_t,UV_t,rho_t,mu,flux_t,V_t);
        data.shape = this.IDC;
        if(this.doSubArea) 
            data.Subspace = v2struct(subArea,accFlux);
        end
    end   
    function dxdt = dx_dt(t,x)
        
        % ignore first row of entries. This is mass in subsystem 
        x        = x(nSpecies+1:end);              
        
        x        = reshape(x,3*M,nSpecies);
        
        y  = x(1:M,:);
        u  = x(M+1:2*M,:);
        v  = x(2*M+1:3*M,:);
        uv = x(M+1:3*M,:);
        
        % convective term; C = v.grad
        C        = sparse(1:m,1:m,u)*Diff.Dy1 + sparse(1:m,1:m,v)*Diff.Dy2;
        
        mu_s     = GetExcessChemPotential(x,t,mu);        
        mu_s(markVinf) = 0;
        
        h_s1    = Diff.Dy1*y - Vext_grad(1:m);
        h_s2    = Diff.Dy2*y - Vext_grad(1+m:end);
                
        h_s1(markVinf)  = 0; 
        h_s2(markVinf)  = 0; 
        
        dydt = -kBT*(Diff.div*uv) - (u.*h_s1  + v.*h_s2);
        duvdt    = - [C*u;C*v] - optsPhys.gamma*uv - Diff.grad*mu_s;
  
        if(doHI)
            error('HI not implemented yet');
%             rho_s    = exp((x-Vext)/kBT);
%             rho_s    = [rho_s;rho_s];
%             gradMu_s = Diff.grad*mu_s;
%             HI_s     = ComputeHI(rho_s,gradMu_s,IntMatrHI);            
%             dxdt     = dxdt + kBT*Diff.div*HI_s + eyes*( h_s.*HI_s );  
        end
        
        
        % need to think about this
        duvdt(Ind.finite,:)     = Ind.normalFinite*uv;       
        
        dydt(markVinf)         = x(markVinf) - x_ic(markVinf);

        dxdt = [dydt;duvdt];
        
        dxdt = [(Int_of_path*GetFlux(x,t))';dxdt(:)];
    end
    function mu_s = GetExcessChemPotential(x,t,mu)
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R) + ...
                                Fex_Meanfield(rho_s,Conv,kBT);
        
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu(iSpecies);
        end

        mu_s = mu_s + x + getVAdd(y1S,y2S,t,optsPhys.V1);
    end   
    function flux = GetFlux(x,t)
        x        = reshape(x,3*M,nSpecies);        
        y  = x(1:M,:);
        uv = x(M+1:3*M,:);
        
        rho_s = exp((y-Vext)/kBT);       
        flux  = -[rho_s;rho_s].*uv;                              
        if(polarShape)
            %then transform to cartesian corrdinates
            flux = GetCartesianFromPolarFlux(flux,ythS);
        end
    end

end