function ShortRangeInteractions_TolmanLength()

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
% 
%% Yukawa-potential
% 
% $$\phi(r) = -\epsilon \frac{e^{-\lambda r}}{4\pi r}$$
% 

    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'TolmanLength'],'ORG');    
    close all;

    %% Parameters               
    epw      = 0.3;
    lambda   = 1;
    lambdaW  = 1;
    R        = inf;
    kBT      = 0.08;    
    
    V2       = struct('V2DV2',@PhiYukawa);
    optsPhys = struct('V2',V2,'HSBulk','CarnahanStarling','kBT',kBT);
    
    N = 200;
    L = 40;
    
    GetCriticalPoint(optsPhys,[],true);    
    disp('Critical values from Sullivan, Wetting Transitions at Fluid-Solid Interfaces');
    disp(['kBT_C = ',num2str(1/11.102)]);
    disp(['rho_C = ',num2str(0.249)]);    
    
    [rhoGas_sat,rhoLiq_sat,mu_sat,p] = BulkSatValues(optsPhys,[0.05;0.5;-0.3],true);    
    
    
    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape     = struct('N',N,'L',L);    
    plotShape = struct('yMin',0,'yMax',30,'N',200);

    IS         = HalfInfSpectralLine(shape);    
    [Pts,Diff] = IS.ComputeAll(plotShape);        
        
    if(R == inf)
        r    = Pts.y;
    else
        r    = Pts.y + R;
    end
    Dr   = Diff.Dy;
    DDr  = Diff.DDy;   
        
    [V,lapV] = Vext(r);
    figure('Name','External Potential');
    subplot(1,2,1);  
    IS.plot(V);
    subplot(1,2,2);  IS.plot(lapV);    
    
    %% Solve ODE
    %
    fig1 = figure('color','white','Position',[0 0 800 600]);
    rho = rhoGas_sat*ones(N,1);
    
    for dmu = -(1e-4):(2e-5):(-2e-5)
        
        [rhoGas_eq,rhoLiq_eq,pLiq,pGas] = ...
            BulkValues(mu_sat+dmu,optsPhys,[rhoGas_sat,rhoLiq_sat,mu_sat],true);            
        rho = fsolve(@ODE,rho);
    
        %hold off;
        IS.plot(rho,'plain'); hold on;
        plot(Pts.y,rhoLiq_eq*ones(size(r)),'b');
        plot(Pts.y,rhoGas_eq*ones(size(r)),'b');
        %title(['dmu = ',num2str(dmu)]);        
        
     %   pause(0.05);
    end
    xlabel('$z$','Interpreter','Latex','fontsize',20);
    ylabel('$\rho$','Interpreter','Latex','fontsize',20);
    
    saveas(fig1,[dirData filesep 'ShortRangePlanarWall.fig']);
   print2eps([dirData filesep 'ShortRange'],fig1);
    
    
    %% Auxiliary functions    
    function y = ODE(rho)  
        
       [muHS,dmuHS,ddmuHS] = mu_HS(rho);       
       drho   = Dr*rho;
       ddrho  = DDr*rho;
       
       if(R == inf)
           LapMu = ddrho.*dmuHS + (drho.^2).*ddmuHS;
       else
           LapMu = 2*drho.*dmuHS./r + ddrho.*dmuHS + (drho.^2).*ddmuHS;
       end
       
       y     = LapMu + rho + lapV - lambda^2*(muHS + V - (mu_sat + dmu));
       
       %Boundary condition
       if(R == inf)
           y(1) = muHS(1) + 2*V(1) - (mu_sat + dmu) - dmuHS(1)*drho(1);
       end
    end
    

    %% Right hand side of ODE    
    function [V,lapV] = Vext(r)
        if(R == inf)
            V    = -epw*exp(-r);
            lapV = -epw*exp(-r);
        else
            V    = -epw/(2*lambdaW^2)*( (R./r - 1./(lambdaW*r)).*exp(-lambdaW*(r-R)) + ...
                                (R./r + 1./(lambdaW*r)).*exp(-lambdaW*(r+R)));
            lapV = (diag(2./r)*Dr + DDr)*V;        
        end
    end

    function [muHS,dmuHS,ddmuHS] = mu_HS(rho)  
         [muHS,~,dmuHS,ddmuHS] = CarnahanStarling(rho,kBT);
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

