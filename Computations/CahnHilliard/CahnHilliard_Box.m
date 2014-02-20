function data = CahnHilliard_Box(optsPhys,optsNum)
%************************************************************************* 
%data = CahnHilliard_Box(optsPhys,optsNum)
%
%            (Eq)  0 = 1/Cn * d(rho*f)/drho -Cn*Lap(rho) - mu + V_ext
%                    = (rho*(1-rho)*(1-2rho))/Cn -Cn*Lap(rho) - mu + V_ext
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%      (BC In/Out) 0 = normal*grad(rho_s);
%
%*************************************************************************   
    if(nargin == 0)
        [optsNum,optsPhys] = Test_CahnHillardBox();
    end
      
    disp(['** ',optsNum.DDFTCode,' **']);
 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
    
    kBT         = optsPhys.kBT; 
    nParticles  = optsPhys.nParticlesS;            
    plotTimes   = optsNum.plotTimes;
    Cn          = optsPhys.Cn;
    Ca          = optsPhys.Ca;
    theta       = optsPhys.theta;
        
    PhysArea    = optsNum.PhysArea;
    N1 = PhysArea.N(1); N2 = PhysArea.N(2);                       
   
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic
    abox                      = Box(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = abox.ComputeAll(optsNum.PlotArea);
        
    %Compute matrix of Integration over subspace    
    subBox        = Box(optsNum.SubArea);    
    IP            = abox.SubShapePts(subBox.Pts);
    Int_SubOnFull = subBox.ComputeIntegrationVector()*IP;
        
    Int_of_path   = subBox.IntFluxThroughDomain(100)*blkdiag(IP,IP);       

%    [Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
%                @SpectralSpectral_Interpolation,Pts,Maps,'normal');        
            
    I       = eye(N1*N2);
    eyes    = [I I];
    %Compute expected curvature
    kappa = 2*cos(theta)/(PhysArea.y1Max - PhysArea.y1Min);
    
    t_preprocess = toc;
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic
    Vext    = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys);
    x_ic    = fsolve(@f_eq,[0;InitialGuess()]);    
    mu      = x_ic(1);
    rho_ic  = x_ic(2:end);
    
    subplot(2,2,1); abox.doPlots(rho_ic,'SC');
    [fE,dfE]    = f(rho_ic);
    subplot(2,2,2); abox.doPlots(-(rho_ic.^2).*dfE/Cn);
    subplot(2,2,3); abox.doPlots(rho_ic.*fE/Cn-mu*rho_ic);
    subplot(2,2,4); abox.doPlots(Diff.Lap*GetExcessChemPotential(rho_ic,0,mu));
    t_eqSol = toc;
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    tic
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('OutputFcn',@PlotFcn,'RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
    [outTimes,X_t] =  ode15s(@drho_dt,plotTimes,[0;rho_ic],opts);        
    t_dynSol = toc;
    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    accFlux   = X_t(:,1);
    X_t       = X_t(:,2:end)';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end
    
    Subspace = v2struct(Path,InterpPath,Int_of_path,Int_SubOnFull,accFlux);
    data     = v2struct(Pts,Diff,Int,Ind,Interp,Conv,X_t,rho_t,Subspace,...
                        mu,flux_t,...
                        t_preprocess,t_eqSol,t_dynSol);                        
    
    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end                    
                    
	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);
    display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
    
    PlotDDFT(v2struct(optsPhys,optsNum,data));    
             
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,x)
                
        rho_s    = x(:,1);
        
        %drhodt   =  
        

        
        %Boundary Conditions: no flux through the walls        
        %flux_dir           = Diff.grad*mu_s;
        %drhodt(Ind.bound)  = Ind.normal*flux_dir;        
%        drhodt(Ind.top)    = rho_s(Ind.top) - rho_ic(Ind.top);%Ind.normalTop*(Diff.grad*rho_s);
%        drhodt(Ind.bottom) = rho_s(Ind.bottom) - rho_ic(Ind.bottom);%Ind.normalBottom*(Diff.grad*rho_s);
             
        dxdt = [drhodt dudt(1:N1*N2) dvdt(N1*N2:end)];
    end

    function y = f_eq(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1);
        rho_s        = x(2:end);
        y            = GetExcessChemPotential(rho_s,0,mu_s); 
        
        %Boundary Conditions
        rhoB         = rho_s(Ind.bound);
        y(Ind.bound) = Cn*(Ind.normal*(Diff.grad*rho_s)) + cos(theta)*(rhoB-1).*rhoB;
        y(Ind.top)   = Ind.normalTop*(Diff.grad*rho_s);
        y(Ind.bottom)= Ind.normalBottom*(Diff.grad*rho_s);
        
        y            = [Int*rho_s - nParticles;y];
    end

    function flux = GetFlux(rho_s,t)
        mu_s  = GetExcessChemPotential(rho_s,t,mu); 
        %flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
        flux  = -Diff.grad*mu_s;                                
    end

    function mu_s = GetExcessChemPotential(rho_s,t,mu_offset)    
        [fE,dfE]    = f(rho_s);
        mu_s        = (rho_s.*dfE + fE)/Cn - Cn*(Diff.Lap*rho_s) + ...
                        - mu_offset + ...
                         getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys);                
    end
    
    function [z,dz] = f(rho)
        z   = 0.5*rho.*(1-rho).^2;
        dz  = 0.5*(1-rho).*(1-3*rho);
    end

    function z = InitialGuess()
        z = 0.5*(1 - tanh( Pts.y2_kv/(2*Cn)));
    end

    %***************************************************************
    %Auxiliary functions:
    %***************************************************************
    function [kBT,nParticles,Cn,Ca,theta] = LoadPhysData_CH(optsPhys)
        kBT         = optsPhys.kBT;
        nParticles  = optsPhys.nParticles;
        Cn          = optsPhys.Cn;
        Ca          = optsPhys.Ca;
        theta       = optsPhys.theta;
    end

    function stat = PlotFcn(t,rho,flag)
        if(~ isscalar(t))
            t = 0;
        end        
        
        stat = 0;
 
        if isempty(flag)
            
            doPlots_IP(Interp,rho);      
            zlim([0 max(rho_ic)]);
            drawnow;
        end
    end

end