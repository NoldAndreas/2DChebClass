function DiffusionHalfInfiniteCapillary()

    disp('** Diffusion HalfInfinite Capillary **');
	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************    
    N1 = 20;
    N2 = 20;    
    L1 = 1;	
    L2 = 1;    
    
    Maps = struct('PhysSpace1',@Comp_to_Phys_1,...
                  'PhysSpace2',@Comp_to_Phys_2,...
                  'CompSpace1',@Phys_to_Comp_1,...
                  'CompSpace2',@Phys_to_Comp_2);
         
    x1PlotMax = Phys_to_Comp_1(3);
    y1Plot    = Comp_to_Phys_1((-1:0.03:x1PlotMax)');
    y2Plot    = Comp_to_Phys_2((-1:0.02:1)');                  
              
    [Pts,Diff,Int,Ind,Interp] = SpectralSpectral(Maps,N1,N2,y1Plot,y2Plot);              
    
    mask1  = (Pts.y1_kv < inf);              
    rho_ic = Vext7(Pts.y1_kv,Pts.y2_kv);
          
    doPlots_IP(Interp,rho_ic);

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
        plot(outTimes(i),Int(mask1)*(rho(mask1)-rho_ic(mask1)),'o');
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
        
        %no flux at origin
        dydt(Ind.left)  = Ind.normalLeft*drho;
        %zero value at infinity
        dydt(Ind.right) = rho(Ind.right);
        %no flux at top and bottom wedge
        dydt(Ind.top & ~Ind.left & ~Ind.right)   = Ind.normalTop(2:end-1,:)*drho;
        dydt(Ind.bottom & ~Ind.left & ~Ind.right)= Ind.normalBottom(2:end-1,:)*drho;                        
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