function theta =IntermediateStokesSlope(z,lambda)
    if(nargin == 0)
        z      = (0.01:0.02:0.2)';
        lambda = 1;        
        boolPlot = true;
        Ca      = 0.05;
    else
        boolPlot = false;
    end        
    
    N          = 150;
    shape      = struct('N',N,'yMin',0,'yMax',pi);    
    plotShape  = struct('yMin',0,'yMax',pi,'N',400);

    SL  = SpectralLine(shape);    
    SL.ComputeAll(plotShape);
    
    y    = SL.Pts.y;  
    IntM = zeros(N);
	vh   = zeros(1,N);
    for i = 2:N
        %Matrix for integral int(v(y),y=y..Y)
        hh          = y(i) - y(i-1);
        vh([i-1,i]) = vh([i-1,i]) + hh/2;
        IntM(i,:)    = vh;
    end
    
    f     = OneOverf(y);    
    
    theta = fsolve(@f_interp,(9*z).^(1/3));
    
    if(boolPlot)
        
        plotOpts = struct('plain',true,'linecolor','b','linewidth',2,'linestyle','-');
        
        figure('color','white','Position',[0 0 800 800]);
        SL.doPlots(IntM*f,plotOpts); hold on; 
        plotOpts.linestyle = ':';
        SL.doPlots(y.^3/9,plotOpts); 

        ylim([0 max(IntM*f)]);
        xlabel('$\theta [rad]$','Interpreter','Latex','fontsize',20);
        ylabel('$G(\theta)$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);
        
        print2eps('theta_G',gcf);  
        saveas(gcf,'theta_G.fig');    
        
        figure('color','white','Position',[0 0 800 800]);
        plot(exp(IntM*f/Ca),y*180/pi,'b','linewidth',2); hold on;
        plot(exp(y.^3/9/Ca),y*180/pi,'b:','linewidth',2);
        xlim([1 max(exp(IntM*f/Ca))])        
        xlabel('$r/\hat L$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta(r)[^\circ]$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);
        
        print2eps('r_theta',gcf);
        saveas(gcf,'r_theta.fig');
        %subplot(1,2,2); plot(theta,z,'ko','markersize',6,'MarkerFaceColor','r');
    end
 
    function f = OneOverf(t)
        f = (lambda*(t.^2-(sin(t)).^2).*(pi-t+sin(t).*cos(t))+((pi-t).^2-(sin(t)).^2).*(t-sin(t).*cos(t)))./ ...
                (2*sin(t).*(lambda^2*(t.^2-(sin(t)).^2)+2*lambda*(t.*(pi-t)+(sin(t)).^2)+(pi-t).^2-(sin(t)).^2));
        f(t==0)  = 0;
        f(t==pi) = 0;
    end

    function h = f_interp(y)
        IP = SL.ComputeInterpolationMatrix(SL.CompSpace(y));
        h  = IP.InterPol*(IntM*f)-z;
    end

end