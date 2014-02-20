function DiffusionClosedBoxFD()

    disp('** DiffusionClosedBoxFD **');
    close all;
    if(length(dbstack) == 1)
        AddPaths();
    end    
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************
    L = 2;
	%Check Polar Spectral/Fourier map in 2D:
    Maps = struct('PhysSpace1',@Comp_to_Phys,...
                  'PhysSpace2',@Comp_to_Phys,...
                  'CompSpace1',@Phys_to_Comp,...
                  'CompSpace2',@Phys_to_Comp);
              
    % no interpolation in FD case
    y1Plot=[];
    y2Plot=[];
    
    N1=20;
    N2=20;
    
    [Pts,Diff,Int,Ind,Interp] = FDFD(Maps,N1,N2,y1Plot,y2Plot);              
        
    [rho_ic] = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
    [rho_ic_IP,Lap_ic_IP] = AnalyticalSolution(Interp.pts1,Interp.pts2,0);
         
    n1     = length(Pts.x1);
    n2     = length(Pts.x2);
 
    doPlots_IP(Interp,rho_ic);
    figure
    doPlots_IP(Interp,Diff.Lap*rho_ic,Lap_ic_IP);

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM  = ones(n1*n2,1);
    mM(Ind.bound) = 0;
    %opts = odeset('OutputFcn',@PlotFcn,'RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.2:5],rho_ic,opts);
    
    %************************************************
    %****************  Postprocess  ****************
    %************************************************
    figure;
    y1M     = reshape(Interp.pts1,Interp.Nplot2,Interp.Nplot1);
    y2M     = reshape(Interp.pts2,Interp.Nplot2,Interp.Nplot1);

    for i=1:length(outTimes)
        rho = Rho_t(i,:)';
           
        z = Interp.InterPol*rho;
        
        subplot(2,1,1)
        mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));        
        title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);       
        zlim([0 max(rho_ic)]);
        
        subplot(2,1,2)
        h = max( rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,outTimes(i)));
        hI = max( z -  AnalyticalSolution(Interp.pts1,Interp.pts2,outTimes(i)));
        hold on;
        plot(outTimes(i),h,'o');
        hold on;        
        plot(outTimes(i),hI,'r');
        
        pause(0.02);
        
    end        
    
    function dydt = Lap(t,rho)
        dydt            = Diff.Lap*rho;
        dydt(Ind.bound) = Ind.normal*(Diff.grad*rho);
    end
                   
    function [y,Lap] = AnalyticalSolution(y1,y2,t)
        y    = cos(pi*y1/L).*cos(pi*y2/L)*exp(-2*t*(pi/L)^2);
        Lap  = -2*(pi/L)^2*cos(pi*y1/L).*cos(pi*y2/L)*exp(-2*t*(pi/L)^2);
    end

    function stat = PlotFcn(t,rho,flag)
        if(~ isscalar(t))
            t = 0;
        end        
        
        stat = 0;
 
        if isempty(flag)
            
            %rho_full = zeros(size(rho_0));
            %rho_full(Ind.outR == 0) = rho;
            rho_full = rho;
            
            doPlots_IP(Interp,rho_full);      
            zlim([0 max(rho_ic)]);
            drawnow;
        end
    end

     %***************************************************************
     %***************************************************************
     %Mapping functions:
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys(xR)
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xR,0,L);
    end
    function xf = Phys_to_Comp(z)
         xf  = -1 + 2*z/L;                     
    end 

    function doPlots_IP(Interp,V,VP)        

        nSpecies=size(V,2);
        if(nSpecies == 1)
            nCol = 1;
        else
            nCol = 2;
        end    
        nRows=ceil(nSpecies/nCol);

        y1 = Interp.pts1;
        y2 = Interp.pts2;   

        if(length(V) == length(Interp.pts1))
            z = V;
        else
            z = real(Interp.InterPol*V);        
        end

        y1M     = reshape(y1,Interp.Nplot2,Interp.Nplot1);
        y2M     = reshape(y2,Interp.Nplot2,Interp.Nplot1);

        xl = [(min(y1)-0.5) (max(y1)+0.5)];
        yl = [(min(y2)-0.5) (max(y2)+0.5)];

        if(nargin == 2)

            if(nSpecies >= 2)

                for iSpecies=1:nSpecies
                    subplot(nRows,nCol,iSpecies);
                    mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));         

                    SetXYLabel();                                
                    title(['Species ' num2str(iSpecies)]);
                end
            else
                mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));         
                SetXYLabel();
            end

        elseif(nargin == 3)

            if(length(VP) == length(Interp.pts1))
                zP = VP;
            else
                zP = real(Interp.InterPol*VP);        
            end    

            subplot(2,1,1)        
            mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));        
            SetXYLabel();

            subplot(2,1,2)        
            mesh(y1M,y2M,reshape(z-zP,Interp.Nplot2,Interp.Nplot1));                           
            SetXYLabel();
        end
        pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);    
        set(gca,'fontsize',15);        

        function SetXYLabel()
            h = xlabel('$y_1$');
            set(h,'Interpreter','Latex'); 
            set(h,'fontsize',20);

            h = ylabel('$y_2$');
            set(h,'Interpreter','Latex'); 
            set(h,'fontsize',20);
        end


    end

end