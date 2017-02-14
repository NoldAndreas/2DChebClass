function DiffusionAdvectionPolarInfinity()
%************************************************
% Diffusion Equation: 
%  drhodt = Lap(rho) - u'*grad(rho);
% Analytical Solution
%  xt = x - ux*t;   y = y - (3 + uy*t);  t0 = -1;
%  f  = exp(-(x^2 + y^2)/(4*(t-t0)))./(t-t0);
%************************************************

    disp('** DiffusionAdvectionPolarInfinity **');
   
    %************************************************
    %****************  Preprocess  ****************
    %************************************************
	N = [24;24];       %collocation points
    L = 3;             %map parameter
    ux  = 0; uy = -1;  %External velocity field:    
    
    IDC                = InfDisc(v2struct(L,N));     
    [Pts,Diff,~,Ind]   = IDC.ComputeAll();    
    
    x1plotMax          = IDC.CompSpace1(10);                      
    IDC.ComputeInterpolationMatrix((0:100)'*x1plotMax/100,(0:0.02:1)',true,true);
                                
    uxV     = ux*ones(IDC.M,1);     
    uyV     = uy*ones(IDC.M,1);
    [ur,ut] = GetPolarFromCartesian(uxV,uyV,Pts.y2_kv);

    rho_ic = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);    
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM           = ones(size(Pts.y1_kv));
    mM(Ind.outR) = 0;    
    opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'MvPattern',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,0:0.2:6,rho_ic,opts);
    
    %****************************************************************
    %***********************  Plot Solution *************************
    %****************************************************************
    Rho_t = Rho_t';                
        
    figure('color','white','Position',[100 100 1100 1100]);    
    for i=1:length(outTimes)        
        rho      = Rho_t(:,i);
        error    = max(abs( rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,outTimes(i))));
        
        IDC.plot(rho,'SC');                
        title(['t = ', num2str(outTimes(i)),' Error: ',num2str(error,'%1.2e')]);
        zlim([0 max(rho_ic)]);
        hold off;                        
    end        
    
    
	%******************************************************
    %*********************** Functions ********************
    %******************************************************
        
    
    function dydt = Lap(t,rho)
        dydt           = Diff.Lap*rho - ([diag(ur);diag(ut)])'*(Diff.grad*rho);
        dydt(Ind.outR) = 0;
    end                   

    function y = AnalyticalSolution(r_o,theta_o,t)
        
        x0 = 0 + ux*t;
        y0 = 3 + uy*t;
        
        [x,y]   = pol2cart(theta_o,r_o);
        y((theta_o == 0) & (r_o == inf)) = 0;
        [theta,r]   = cart2pol(x-x0,y-y0);
        
        t0 = -1;
        y  = exp(-r.^2./(4*(t-t0)))./(t-t0);
    end

end