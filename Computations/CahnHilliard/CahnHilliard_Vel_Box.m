function data = CahnHilliard_Vel_Box(optsPhys,optsNum)
%NOT FULLY ADAPTED YET
%************************************************************************* 
%data = CahnHilliard_Vel_Box(optsPhys,optsNum)
%
%          (Eq)  0 = 1/Cn * d(rho*f)/drho -Cn*Lap(rho) - mu - V_ext
%                  = (rho*(1-rho)*(1-2rho))/Cn -Cn*Lap(rho) - mu - V_ext
%          (Dyn 1) drho/dt + div(u*rho) = 0
%          (Dyn 2) rho*(du/dt + (u*grad)u) = - grad(p) + ...
%                                          + eta*Lap(u) + deta/dxj*dui/dxj
%                                          + zeta*grad(div(u))+
%                                          + div(u)*grad(zeta)
%                                          + Cn/Cak*rho*grad(Lap(rho))
%            where p = rho^2 * f'(rho)/(Cn*Cak) 
%   Here, we assume that 
%           (a)eta  = etaRho*rho and 
%           (a)zeta = zetaRho*rho  
%*************************************************************************   
    if(nargin == 0)
        %Numerical Parameters    
        Phys_Area = struct('y1Min',0,'y1Max',10,'N',[22,22],...
                           'y2Min',-10,'y2Max',10,'L2',4);

        Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                           'y2Min',-10,'y2Max',10,'L2',3,'N2',100);

        %Sub_Area  = Phys_Area;
        Sub_Area = struct('y1Min',5,'y1Max',10,'N1',20,...
                          'y2Min',0,'y2Max',5,'N2',20);    

        optsNum = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                         'DDFTCode','CahnHilliard_Vel_Box',...
                         'plotTimes',0:0.2:20);                            

        optsPhys = struct('V1DV1','Vext_Cart_TestCahnHillard_1',...
                          'epsilon_w',0,'epsilon_w_end',0.0,...
                          'theta',pi/4,'Cn',0.6,...
                          'Ca',0.01,'Cak',0.6,'etaRho',1,'zetaRho',1,...
                          'kBT',0.7,...
                          'nParticlesS',100);
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
    Cak         = 6*optsPhys.Ca;
    theta       = optsPhys.theta;
    etaRho      = optsPhys.etaRho;
    zetaRho     = optsPhys.zetaRho;    
        
    PhysArea    = optsNum.PhysArea;
    N1 = PhysArea.N(1); N2 = PhysArea.N(2);                       
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic
    abox                      = Box(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = abox.ComputeAll(optsNum.PlotArea);    

%    [Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
%                @SpectralSpectral_Interpolation,Pts,Maps,'normal');        
            
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
    
    subplot(2,2,1); abox.doPlots(rho_ic);
    [fE,dfE]    = f(rho_ic);
    subplot(2,2,2);  abox.doPlots(-(rho_ic.^2).*dfE/Cn);
    subplot(2,2,3);  abox.doPlots(rho_ic.*fE/Cn-mu*rho_ic);
    subplot(2,2,4);  abox.doPlots(Diff.Lap*GetExcessChemPotential(rho_ic,0,mu));
    
    t_eqSol = toc;
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    tic
    mM             = ones(N1*N2,1);
    mM(Ind.bound)  = 0;
    opts           = odeset('RelTol',10^-8,'AbsTol',10^-8,...
                            'Mass',diag(repmat(mM,3,1)));
    x_ic           = [rho_ic ; zeros(2*N1*N2,1)];
    [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);        
    
    t_dynSol = toc;
    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    X_t       = X_t';
    rho_t     = X_t(1:N1*N2,:);
    
    flux_t    = X_t(N1*N2+1:end,:);
    
    data     = v2struct(Pts,Diff,Int,Ind,Interp,X_t,rho_t,...
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
        
        x        = reshape(x,N1*N2,3);
        rho      = x(:,1);
        u        = x(:,2);
        v        = x(:,3);
        uv       = x(:,2:3);
        uv       = uv(:);
        rho2     = repmat(rho,2,1); 

        drhodt   =  - Diff.div*(uv.*rho2);
        
        %Convective term:        
        C        = [diag(u) diag(v)]*Diff.grad;
        CC       = blkdiag(C,C);
        
        
        duvdt    =  - CC*uv +... 
                    - (1 - 6*rho2 + 6*rho2.^2).*(Diff.grad*rho)/(Cn*Cak) + ...
                    + etaRho*(Diff.LapVec*uv) + ...
                    + zetaRho*(Diff.gradDiv*uv) + ...
                    + (Diff.gradLap*rho)*Cn/Cak;                            
        
        %Pressure Term = - grad(p)/rho 
        %              = - 1/rho*d(rho^2 *df )/drho *grad(rho)/(Cn*Cak)
        %              = - (1 - 6rho + 6rho^2)*grad(rho)/(Cn*Cak)
                
        rhoL         = rho(Ind.left);
        rhoR         = rho(Ind.right);
        %Boundary Conditions for velocities
        duvdt([Ind.left;Ind.left])    = uv([Ind.left;Ind.left]);
        duvdt([Ind.right;Ind.right])  = uv([Ind.right;Ind.right]);
        [uvBottom,uvTop] = FlowAtEntries(t);
        duvdt([Ind.top;Ind.top])      = uvTop - uv([Ind.top;Ind.top]);
        duvdt([Ind.bottom;Ind.bottom])= uvBottom - uv([Ind.bottom;Ind.bottom]);        
        %duvdt([Ind.top;Ind.top])      = uv([Ind.top;Ind.top]);
        %duvdt([Ind.bottom;Ind.bottom])= uv([Ind.bottom;Ind.bottom]);        
        
        %Boundary Conditions for density
        drhodt(Ind.left)  = Cn*(Ind.normalLeft*(Diff.grad*rho)) + cos(theta)*(rhoL-1).*rhoL;
        drhodt(Ind.right) = Cn*(Ind.normalRight*(Diff.grad*rho)) + cos(theta)*(rhoR-1).*rhoR;
        drhodt(Ind.top)   = Ind.normalTop*(Diff.grad*rho);
        drhodt(Ind.bottom)= Ind.normalBottom*(Diff.grad*rho);

        %flux_dir           = Diff.grad*mu_s;
        %drhodt(Ind.bound)  = Ind.normal*flux_dir;        
%        drhodt(Ind.top)    = rho_s(Ind.top) - rho_ic(Ind.top);%Ind.normalTop*(Diff.grad*rho_s);
%        drhodt(Ind.bottom) = rho_s(Ind.bottom) - rho_ic(Ind.bottom);%Ind.normalBottom*(Diff.grad*rho_s);
             
        dxdt = [drhodt ; duvdt];
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
    
    function [uvBottom,uvTop] = FlowAtEntries(t)
        Z = zeros(N1*N2,1);
        tau   = 2;
        uMaxt = -0.1*(1- exp(-(t/tau)^2));
        h     = uMaxt*(1 - (2*Pts.y1_kv/PhysArea.y1Max - 1).^2);
        
        uvBottom = [Z(Ind.bottom) ; h(Ind.bottom)];
        uvTop    = [Z(Ind.top) ; h(Ind.top)];
    end

    %***************************************************************
    %Mapping functions:
    %***************************************************************
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys1(x)
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(x,PhysArea.y1Min,PhysArea.y1Max);
    end
    function xf = Phys_to_Comp1(z)
        xf = InvLinearMap(z,PhysArea.y1Min,PhysArea.y1Max);        
    end

    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys2(x)    
        [z,dz,dx,ddx,dddx,ddddx] = QuadMapAB(x,PhysArea.L2,PhysArea.y2Min,PhysArea.y2Max);
    end

    function xf = Phys_to_Comp2(z)           
        xf = InvQuadMapAB(z,PhysArea.L2,PhysArea.y2Min,PhysArea.y2Max);        
    end



end