function DiffusionHalfSpace()

    disp('** Diffusion Half Space **');
	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************    
    N1 = 30;
    N2 = 30;    
    L1 = 2;	L2 = 2;    
    PlotArea = struct('y1Min',-5,'y1Max',5,'N1',100,...
                      'y2Min',0,'y2Max',10,'N2',100);
    
    HS = HalfSpace(v2struct(L1,L2),[N1;N2]);
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(PlotArea);
            
    rho_ic = Vext15(Pts.y1_kv,Pts.y2_kv);
          
    doPlots_SC(Interp,Pts,rho_ic);

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************        
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.1:5],rho_ic,opts);

    for i=1:length(outTimes)
        
        rho = Rho_t(i,:)';           
        
        subplot(2,2,1)
        doPlots_IP(Interp,rho);        
        title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);       
        zlim([0 max(rho_ic)]);
        
        subplot(2,2,2)
        hold on;
        plot(outTimes(i),Int*(rho-rho_ic),'o');
        title('Gain of Mass');
        
        subplot(2,2,3);
        fl   = -Diff.grad*rho;                
        fl_x = fl(1:N1*N2);
        fl_y = fl(N1*N2+1:end);
        quiver(Phys_to_Comp_1(Pts.y1_kv),...
            Phys_to_Comp_2(Pts.y2_kv),fl_x,fl_y);        
        title(['Flux, max: ' , num2str(max(fl_x.^2 + fl_y.^2))]);        
        
        pause(0.05);
        
    end
    
    function dydt = Lap(t,rho)
        dydt             = Diff.Lap*rho;
        
        drho             = Diff.grad*rho;
        
        dydt(Ind.bound)  = rho(Ind.bound);
        dydt(Ind.bottom) = Ind.normalBottom*drho;
        
    end                   

    %***************************************************************
    %   Mapping functions:
    %***************************************************************         
    function [z,dz,dx,ddx,dddx] = Comp_to_Phys_1(xR)
        [z,dz,dx,ddx,dddx] = PosRay(xR,L1);
    end
    function [z,dz,dx,ddx,dddx] = Comp_to_Phys_2(xT)    
        [z,dz,dx,ddx,dddx] = LinearMap(xT,-L2,L2);
    end
    function xf = Phys_to_Comp_1(z)
         xf  = InvPosRay(z,L1);
    end
    function xf = Phys_to_Comp_2(z)           
        xf = z/L2;
    end

end