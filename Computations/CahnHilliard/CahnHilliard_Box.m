function CahnHilliard_Box()
    
    h = 5;
 
	PhysArea = struct('y1Min',-h,'y1Max',h,'N',[40,40],...
                       'y2Min',0,'y2Max',h);
                   
	PlotArea = PhysArea;
    PlotArea.N1 = 100;  PlotArea.N2 = 100;            
                              
    optsPhys = struct('theta',pi/4,'Cn',1,...                      
                      'rho_m','0',...
                      'nParticlesS',0);          
     
    %%  Initialization
    close all;
        
    nParticles  = optsPhys.nParticlesS;            
    Cn          = optsPhys.Cn;    
    theta       = optsPhys.theta;
    rho_m       = optsPhys.rho_m;
        
    abox               = Box(PhysArea);
    [Pts,Diff,Int,Ind] = abox.ComputeAll(PlotArea);
        
    %Compute expected curvature
    kappa = 2*cos(theta)/(PhysArea.y2Max - PhysArea.y2Min);
            
    %% Compute Equilibrium    
    x_ic    = fsolve(@f_eq,[0;InitialGuess()]);    
    mu      = x_ic(1);
    rho     = x_ic(2:end);
    
    
    [dW,W]         = DoublewellPotential(rho,Cn);    
    p              = -W+mu*(rho + rho_m);
    p              = p-mean(p(Ind.left));
    surfaceTension = 4/3;
    deltaP         = mean(p(Ind.right)-p(Ind.left));
    
    %% Postprocess
    disp(['Pressure difference from curvature: ',num2str(kappa*surfaceTension)]);
    disp(['Pressure difference - measured: ',num2str(deltaP)]);    
    disp(['Error: ',num2str(kappa*surfaceTension - deltaP)]);   
    
    figure('color','white','Position',[0 0 800 600]);    
    abox.plot(rho); hold on;    
    abox.plot(p);
    
    print2eps(['Capillary_',num2str(h)],gcf);        
	saveas(gcf,['Capillary_',num2str(h),'.fig']);    

	figure('color','white','Position',[0 0 600 600]);    
    abox.doPlotFLine([PhysArea.y1Min PhysArea.y1Max],h/2*[1 1],rho,struct('color','b')); hold on;    
    abox.doPlotFLine([PhysArea.y1Min PhysArea.y1Max],h/2*[1 1],p,struct('color','r')); hold on;        
    xlabel('$y_1$','Interpreter','Latex','fontsize',20);
    ylabel('$\varrho,p-p_V$','Interpreter','Latex','fontsize',20);
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',15);
    
    print2eps(['Capillary_CrossSection_',num2str(h)],gcf);        
	saveas(gcf,['Capillary_CrossSection_',num2str(h),'.fig']);    
    
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function y = f_eq(x)                    
        mu_s         = x(1);
        rho_s        = x(2:end);
        y            = GetExcessChemPotential(rho_s,mu_s); 
        
        %Boundary Conditions
        rhoB         = rho_s(Ind.bound);
        y(Ind.bound) = Cn*(Ind.normal*(Diff.grad*rho_s)) + cos(theta)*(rhoB.^2-1);
        y(Ind.left)  = Ind.normalLeft*(Diff.grad*rho_s);
        y(Ind.right)= Ind.normalRight*(Diff.grad*rho_s);
        
        y            = [Int*rho_s - nParticles;y];
    end

    function mu_s = GetExcessChemPotential(rho_s,mu_offset)    
        dW     = DoublewellPotential(rho_s,Cn);
        mu_s   = dW - Cn*(Diff.Lap*rho_s) - mu_offset;                                   
    end
    
    function z = InitialGuess()
        z = tanh(Pts.y1_kv/(2*Cn));
    end

    
end