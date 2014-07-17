function ComputeDynamics(this)
%***************************************************************    
    %***************************************************************
    %   Dynamics:
    %***************************************************************
    %***************************************************************    	
    %****************************************
	I     = eye(N1*N2);   
    eyes  = repmat(I,1,2);    

    
    disp('Compute dynamic 2D evolution ..');         
    mark             = true(N1*N2,1);   
    optsNum.maxComp_y2 = inf;
    opts             = PhysArea;
    opts.optsPhys    = optsPhys;     
    opts.optsNum     = optsNum;     
    opts.optsNum     = rmfield(opts.optsNum,'PhysArea');     
    opts.Comments    = this.configName;
    [sol,rec,paramsDyn]   = DataStorage('DynamicSolutions',@ComputeDynamics,opts,x_ic);                     
    plotTimes         = sol.outTimes;                

    %***************************************************************
    %*******************Postprocess Dynamics************************
    %***************************************************************
    
    XRho_t(mark,:)  = (sol.X_t(:,1:N1*N2))';
	rho_t           = exp((XRho_t-Vext*ones(1,length(plotTimes)))/kBT);
    if(optsPhys.Inertial)         
         flux_t      = [rho_t;rho_t].*(sol.X_t(:,1+N1*N2:end))';
    else
        flux_t          = zeros(2*N1*N2,length(plotTimes));    
        for i = 1:length(plotTimes)
            flux_t(:,i) = GetFlux(XRho_t(:,i),plotTimes(i));
        end
    end            
        
    shapeSub      = struct('y1Min',-2,'y1Max',2,'y2Min',R,'y2Max',3,'N',[40,40]);
    subBox        = Box(shapeSub);
    IP            = HS.SubShapePtsCart(subBox.Pts);
    Int_SubOnFull = subBox.ComputeIntegrationVector()*IP;
        
    Int_of_path   = subBox.IntFluxThroughDomain(100)*blkdiag(IP,IP);
    
    accFlux = zeros(length(plotTimes),1);    
    for i = 2:length(plotTimes)        
        accFlux(i) = accFlux(i-1)+0.5*(plotTimes(i)-plotTimes(i-1))*...
            (Int_of_path*(flux_t(:,i-1)+flux_t(:,i)));
    end    
    
    Subspace      = v2struct(Int_of_path,Int_SubOnFull,accFlux);%InterpPath,Path
    Subspace.subArea = subBox;
    %***************************************************************
    %***********************Plotting********************************
    %***************************************************************
    optsPlot.nContours = [0.1,0.3,0.5,0.7];                
    data       = v2struct(rho_t,flux_t,optsPlot,Subspace);
    data.shape = HS;
    optsNum.plotTimes = plotTimes; 

    plotData = v2struct(optsPhys,optsNum,data);
         
    if(PersonalUserOutput)        
        PlotDDFT(plotData,[dirData filesep 'DynamicSolutions' filesep paramsDyn.Filename(1:end-4)]);            
        PlotDDFT(plotData);
    end

    diary
    dirData = orgdirData;
    %***************************************************************
    %***********************End of main code************************
    %***************************************************************
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************      
    function sol = ComputeDynamics(params,x_ic)
        plotTimes  = (0:optsNum.plotTimes_Interval:optsNum.plotTimes_T);
        if(optsPhys.Inertial)
            mM            = ones(3*N1*N2,1);    
            mM([Ind.top;false(N1*N2,1);Ind.bottom]) = 0;            
            mM            = mM([mark;mark;mark]);
            opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
            [sol.outTimes,sol.X_t] =  ode15s(@dx_dtInertial,plotTimes,...
                                [x_ic(mark);zeros(2*sum(mark),1)],opts);
        else
            mM             = ones(N1*N2,1);    
            mM(Ind.bottom|Ind.top) = 0;
            mM             = mM(mark);
            opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
            [sol.outTimes,sol.X_t] =  ode15s(@dx_dt,plotTimes,x_ic(mark),opts);                
        end
    end    
    function mu_s = GetExcessChemPotentialFull(x)                
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - (mu_sat + Dmu);
        end
        mu_s = mu_s + x + Conv*rho_s + VAdd;           
    end   
    function dxdt = dx_dt(t,x)
       
        VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        
        x(mark)  = x;
        x(~mark) = x_ic(~mark);
        %x = x';
        
        mu_s     = GetExcessChemPotentialFull(x);        
        mu_s((Pts.y2_kv==inf),:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        dxdt     = dxdt/optsPhys.gamma;
                
        %Boundary Conditions at infinity
        flux_dir           = Diff.grad*mu_s;
        dxdt(Ind.bottom,:) = Ind.normalBottom*flux_dir;
        dxdt(Ind.top,:)    = x(Ind.top)-x_ic(Ind.top);

        dxdt = dxdt(:);      
        dxdt = dxdt(mark);
    end
    function dxdt = dx_dtInertial(t,x)         
        VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        
        x        = reshape(x,sum(mark),3);
        x(mark,:)  = x;
        x(~mark,1) = x_ic(~mark);  %outside computational domain, set density to intial and set velocity to zero
        
        y        = x(:,1);
        u        = x(:,2);
        v        = x(:,3);
        uv       = x(:,2:3);
        uv       = uv(:);
        rho      = exp((y-Vext)/kBT);        
        
        %Convective term:        
        %C        = [diag(u) diag(v)]*Diff.grad;
        %C        = sparse(diag(u))*Diff.Dy1 + sparse(diag(v))*Diff.Dy2;
        m        = N1*N2;
        C        = sparse(1:m,1:m,u)*Diff.Dy1 + sparse(1:m,1:m,v)*Diff.Dy2;
        %CC      = blkdiag(C,C);
        mu_s     = GetExcessChemPotentialFull(y);
        mu_s(~mark) = 0;
        mu_s((Pts.y2_kv==inf),:) = 0; %necessary?
        
        %Equations       
        h_s1    = Diff.Dy1*y - Vext_grad(1:m);
        h_s2    = Diff.Dy2*y - Vext_grad(1+m:end);
                
        h_s1(Pts.y2_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y2_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y1_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y1_kv==-inf) = 0; %here, we have assumed that grad(mu) converges fast enough
        
        %h_s      = Diff.grad*y - Vext_grad;
        %h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        %h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        %h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
        
        dydt     = - kBT*(Diff.div*uv) - (u.*h_s1  + v.*h_s2);        
        %dydt     = - kBT*(Diff.div*uv) - [diag(u) diag(v)]*h_s;        
        %duvdt    = - CC*uv - optsPhys.gamma*uv - Diff.grad*mu_s;
        duvdt    = - [C*u;C*v] - optsPhys.gamma*uv - Diff.grad*mu_s;
               
        %Boundary Conditions for velocities
        dydt(Ind.top)                       = y(Ind.top)- x_ic(Ind.top); %Density at y_2 = infinity        
        duvdt([false(N1*N2,1);Ind.bottom])  = uv([false(N1*N2,1);Ind.bottom]); %Normal Velocity at wall       
        dxdt = [dydt(mark); duvdt([mark;mark])];       
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        VAdd  = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        mu_s  = GetExcessChemPotentialFull(x); 
        mu_s((Pts.y2_kv==inf),:) = 0;
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s)/optsPhys.gamma;                                
    end

end