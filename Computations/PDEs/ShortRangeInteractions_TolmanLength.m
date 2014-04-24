function ShortRangeInteractions_TolmanLength()

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
% 
%% Yukawa-potential
% 
% $$\phi(r) = -\epsilon \frac{e^{-\lambda r}}{4\pi r}$$
% 

    %% Parameters               
    epw      = 5;
    lambda   = 0.3;
    lambdaW  = 0.2;
    R        = 10;
    kBT      = 0.75;
    dmu      = -0.1;
    
    V2       = struct('V2DV2',@PhiYukawa);
    optsPhys = struct('V2',V2,'HSBulk','MuCarnahanStarling','kBT',kBT);
    
    N = 200;
    L = 20;
    
    GetCriticalPoint(optsPhys,[],true);    
    [rhoGas_sat,rhoLiq_sat,mu_sat,p] = BulkSatValues(optsPhys,[0.05;0.5;-3],true);    
    
    [rhoGas_eq,rhoLiq_eq,pLiq,pGas] = BulkValues(mu_sat-0.01,optsPhys,[0.05;0.5;-0.3],true);
    
    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape     = struct('N',N,'L',L);    
    plotShape = struct('yMin',0,'yMax',30,'N',200);

    IS         = HalfInfSpectralLine(shape);    
    [Pts,Diff] = IS.ComputeAll(plotShape);        
        
    r    = Pts.y + R;
    Dr   = Diff.Dy;
    DDr  = Diff.DDy;   
        
    [V,lapV] = Vext(r);
    figure('Name','External Potential');
    subplot(1,2,1);  IS.doPlots(V);
    subplot(1,2,2);  IS.doPlots(lapV);    
    
    %% Solve ODE
    %
    figure;    
    rho = rhoGas_sat*ones(N,1);
    
    for dmu = -0.001:0.0001:0.000001
        rho = fsolve(@ODE,rho);
    
        hold off;
        IS.doPlots(rho); hold on;
        plot(r-R,rhoLiq_eq*ones(size(r)),'b');
        plot(r-R,rhoGas_eq*ones(size(r)),'b');
        
        pause(0.05);
    end
    
    %% Auxiliary functions    
    function y = ODE(rho)  
        
       [muHS,dmuHS,ddmuHS] = mu_HS(rho);       
       drho   = Dr*rho;
       ddrho  = DDr*rho;
       
       LapMu = 2*drho.*dmuHS./r + ddrho.*dmuHS + (drho.^2).*ddmuHS;
       
       y     = LapMu + rho + lapV - lambda^2*(muHS + V - (mu_sat + dmu));
       
    end
    

    %% Right hand side of ODE
    function [V,lapV] = Vext(r)
        V    = -epw/(2*lambdaW^2)*( (R./r - 1./(lambdaW*r)).*exp(-lambdaW*(r-R)) + ...
                                (R./r + 1./(lambdaW*r)).*exp(-lambdaW*(r+R)));
        lapV = (diag(2./r)*Dr + DDr)*V;
    end

    function [muHS,dmuHS,ddmuHS] = mu_HS(rho)  
         [muHS,~,dmuHS,ddmuHS] = MuCarnahanStarling(rho,kBT);
         muHS   = muHS   + kBT*(log(rho));
         dmuHS  = dmuHS  + kBT./rho;
         ddmuHS = ddmuHS - kBT./(rho.^2);
    end

    function [z,dzdr_r,alpha] = PhiYukawa(r) 
        % only third output is used               
        %%
        % 
        % $$\alpha = 2\pi \int_0^\infty r^2 \phi(r) dr$$
        %         
        alpha  = -1/(2*lambda^2);
        z      = NaN;
        dzdr_r = NaN;
    end
end

