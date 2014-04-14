function CahnHilliard_InfCapillary()

%************************************************************************* 
%
%  Contact angle always theta = pi/2 (meaning no pressure jump across the interface)
%
%            (Eq)  0 = dW/drho - Cn*Lap(rho) - mu 
%        with W(rho) = (1-rho^2)^2/(2*Cn)
%
% with   (BC Wall) 0 = Cn*normal*grad(rho_s) + cos(theta)*(rho-1)*rho
%                    = Cn*normal*grad(rho_s)
%  (BC +/- infinity) 0 = normal*grad(rho_s);
%
% Dynamic equation (stationary, Reynolds number Re = 0)
%
%%
% 
% $$ 0 = \nabla \cdot ( (\rho + \rho_m){\bf u} )$$
%
% $$ 0 = - (\rho + \rho_m) \nabla \mu + Ca  \nabla \left\{  (\zeta + \eta)(\nabla \cdot {\bf u}){\bf I} + \eta \nabla {\bf u}  \right\} $$
% 
%  
%*************************************************************************   

    %% Parameters
    
    % Numerical Parameters    
    PhysArea  = struct('y2Min',0,'y2Max',1,'N',[40,40],'L1',1);
    PlotArea  = struct('y2Min',0,'y2Max',1,'N1',80,...
                       'y1Min',-1.5,'y1Max',1,'N2',80);                                       
    M = PhysArea.N(1)*PhysArea.N(2);   
             
    % Physical Parameters
    mu_s  = 0;   % Chemical potential (saturation = 0)
    Cn    = 0.4;   % Cahn number (=1,microscopic scaling)
    theta = pi/2; %Equilibrium contact angle
    Ca    = 0.3;   % Capillary number
    rho_m = 2; % = (rhoL+rhoV)/(rhoL-rhoV)
    Uwall = 0.0001;
    Re    = 1;
    
    etaRho  = 1;
    zetaRho = 10;
    
    plotTimes = 0:0.1:20;
                  
    %% Initialization
    
	IC                        = InfCapillary(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = IC.ComputeAll(PlotArea);    
    
    %% Equilibrium
    rho_ic    = fsolve(@f_eq,InitialGuess());    
    IC.doPlots(rho_ic,'SC');
    
    
    %% Dynamics   
    mM             = ones(M,1);
    mM(Ind.top | Ind.bottom)  = 0;
    opts           = odeset('RelTol',10^-8,'AbsTol',10^-8,...
                            'Mass',diag(repmat(mM,3,1)));
    x_ic           = [rho_ic ; zeros(2*M,1)];
    [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);        

    X_t       = X_t';
	rho_t     = X_t(1:M,:);    
    flux_t    = X_t(M+1:end,:);
    %x_iguess = [rho_ic ; Uwall*(Pts.y2_kv/max(Pts.y2_kv)) ; zeros(M,1)];   
    %x        = fsolve(@f_dyn,x_iguess);    
    for i = 1:length(outTimes)
        subplot(1,2,1);
        IC.doPlots(rho_t(:,i),'contour'); hold on;
        IC.doPlotsFlux(flux_t(:,i)); title(['t=',num2str(outTimes(i)),' , max = ',num2str(max(abs(flux_t(:,i))))]); 
        
        subplot(1,2,2);
        IC.doPlots(rho_t(:,i),'SC');
        
        pause(0.1);
    end

    %% Functions    
	function y = f_eq(rho_s)
        y            = GetExcessChemPotential(rho_s,0,mu_s); 
        
        % Boundary Conditions
        rhoB          = rho_s(Ind.bottom);
        y(Ind.bottom) = Cn*(Ind.normalBottom*(Diff.grad*rho_s)) + cos(theta)*(rhoB-1).*rhoB;
        
        rhoB        = rho_s(Ind.top);
        y(Ind.top)  = Cn*(Ind.normalTop*(Diff.grad*rho_s)) + cos(theta)*(rhoB-1).*rhoB;        
              
    end

	function y = f_dyn(x)        
        x        = reshape(x,M,3);
        rho      = x(:,1);
        %u        = x(:,2);
        %v        = x(:,3);
        uv       = x(:,2:3);
        uv       = uv(:);
        rho2     = repmat(rho,2,1); 

        %Continuity
        cont     =  - Diff.div*(uv.*(rho_m+rho2));
        
        % Momentum Equation
        momChemPot  = -(rho2+rho_m).*( -2*(1-3*rho2.^2).*(Diff.grad*rho)/Cn  ...
                                       - Cn*Diff.gradLap*rho);                 
        mom         =  momChemPot + ...
                        + Ca*(etaRho*(Diff.LapVec*uv) +...
                              zetaRho*(Diff.gradDiv*uv));  
        
        %Boundary Conditions for velocities                                  
        mom([Ind.top;Ind.top])      = uv([Ind.top;Ind.top]) - ...
                              [Uwall*ones(sum(Ind.top),1);zeros(sum(Ind.top),1)];
        mom([Ind.bottom;Ind.bottom])= uv([Ind.bottom;Ind.bottom]);        
        
        %Boundary Conditions for density       
        rhoB     = rho(Ind.bottom);     rhoT     = rho(Ind.top);      
        cont(Ind.bottom) = Cn*(Ind.normalBottom*(Diff.grad*rho)) + cos(theta)*(rhoB-1).*rhoB;
        cont(Ind.top)    = Cn*(Ind.normalTop*(Diff.grad*rho)) + cos(theta)*(rhoT-1).*rhoT;        
        
        %drhodt(Ind.top)   = Ind.normalTop*(Diff.grad*rho);
        %drhodt(Ind.bottom)= Ind.normalBottom*(Diff.grad*rho);

        %flux_dir           = Diff.grad*mu_s;
        %drhodt(Ind.bound)  = Ind.normal*flux_dir;        
%        drhodt(Ind.top)    = rho_s(Ind.top) - rho_ic(Ind.top);%Ind.normalTop*(Diff.grad*rho_s);
%        drhodt(Ind.bottom) = rho_s(Ind.bottom) - rho_ic(Ind.bottom);%Ind.normalBottom*(Diff.grad*rho_s);
             
        y = [cont ; mom];                 
    end

    function dxdt = dx_dt(t,x)
        
        x        = reshape(x,M,3);
        rho      = x(:,1);
        uv       = x(:,2:3);
        uv       = uv(:);
        rho2     = repmat(rho,2,1); 

        %Continuity
        drhodt   =  - Diff.div*(uv.*(rho2+rho_m));
        
        % Momentum balance (no nonlinear term)
        duvdt     = -( -2*(1-3*rho2.^2).*(Diff.grad*rho)/Cn  ...
                                       - Cn*Diff.gradLap*rho)/(Ca*Re) ...
                        + 1/Re*(etaRho*(Diff.LapVec*uv) +...
                              zetaRho*(Diff.gradDiv*uv))./(rho2+rho_m);
        
        %Nonlinear term
        %C        = [diag(u) diag(v)]*Diff.grad;
        %CC       = blkdiag(C,C);
        %- CC*uv                
                
        %Boundary Conditions for velocities
        duvdt([Ind.top;Ind.top])      = uv([Ind.top;Ind.top]);
        duvdt([Ind.bottom;Ind.bottom])= uv([Ind.bottom;Ind.bottom]);
        
        %Boundary Conditions for density
        rhoB     = rho(Ind.bottom);     rhoT     = rho(Ind.top);      
        drhodt(Ind.bottom) = Cn*(Ind.normalBottom*(Diff.grad*rho)) + cos(theta)*(rhoB-1).*rhoB;
        thetaUp = UpperTheta(t);
        drhodt(Ind.top)    = Cn*(Ind.normalTop*(Diff.grad*rho)) + cos(thetaUp)*(rhoT-1).*rhoT;        
        

        %flux_dir           = Diff.grad*mu_s;
        %drhodt(Ind.bound)  = Ind.normal*flux_dir;        
%        drhodt(Ind.top)    = rho_s(Ind.top) - rho_ic(Ind.top);%Ind.normalTop*(Diff.grad*rho_s);
%        drhodt(Ind.bottom) = rho_s(Ind.bottom) - rho_ic(Ind.bottom);%Ind.normalBottom*(Diff.grad*rho_s);
             
        dxdt = [drhodt ; duvdt];
    end


    function mu_s = GetExcessChemPotential(rho_s,t,mu_offset)    
        [WE,dWE]    = W(rho_s);
        mu_s        = dWE - Cn*(Diff.Lap*rho_s) - mu_offset;                           
    end

    function [z,dz] = W(rho)        
        z   = (1-rho.^2).^2/(2*Cn);
        dz  = -2*rho.*(1-rho.^2)/Cn;
    end

    function z = InitialGuess()
        z = tanh( Pts.y1_kv/(2*Cn));
        %z = ones(M,1);
    end

    function PlotX(x)
        x        = reshape(x,M,3);
        rhos      = x(:,1);
        us        = x(:,2);
        vs        = x(:,3);
        
        figure;
        subplot(1,2,1);     IC.doPlots(rhos);
        subplot(1,2,2);     IC.doPlotsFlux([us;vs]);
    end

    function th = UpperTheta(t)
        th = max(90 + 10*(t/10),100)*pi/180;
    end

end