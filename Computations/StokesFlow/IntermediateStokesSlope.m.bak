function theta = IntermediateStokesSlope(z,lambda)
    if(nargin == 0)
        z      = (0.01:0.02:0.2)';
        lambda = 0;        
        boolPlot = true;
        Ca      = 0.05;
    else
        boolPlot = false;
    end        
    
    N          = 200;
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
    
    f     = OneOverf(y,lambda);    
    
    theta = fsolve(@f_interp,(9*z).^(1/3));
    
    if(boolPlot)
        
        plotOpts = struct('plain',true,'linecolor','b','linewidth',2,'linestyle','-');
        f1 =  figure('color','white','Position',[0 0 800 800]);
        plot(y,GHR(y),'b','linewidth',2); hold on;
        plotOpts.linecolor = 'r';
        SL.doPlots(IntM*OneOverf(y,1),plotOpts); hold on; 
        plotOpts.linecolor = 'm';
        SL.doPlots(IntM*OneOverf(y,10),plotOpts); hold on; 
        plotOpts.linecolor = 'k';
        plotOpts.linestyle = ':';
        SL.doPlots(y.^3/9,plotOpts);
        ylim([0 max(IntM*f)]);
        set(gca,'xtick',[0,pi/2,pi],...
                'xticklabel',{'0','90','180'});        
        xlabel('$\theta [deg]$','Interpreter','Latex','fontsize',20);
        ylabel('$G(\theta)$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);                
        
        plotOpts = struct('plain',true,'linecolor','b','linewidth',2,'linestyle','-');        
        f2 =  figure('color','white','Position',[0 0 800 800]);
        plot(y,GHR(y),'b','linewidth',2); hold on;
        plotOpts.linecolor = 'r';
        SL.doPlots(IntM*OneOverf(y,1),plotOpts); hold on; 
        plotOpts.linecolor = 'm';
        SL.doPlots(IntM*OneOverf(y,10),plotOpts); hold on; 
        plotOpts.linecolor = 'k';
        plotOpts.linestyle = ':';
        SL.doPlots(y.^3/9,plotOpts); 

        ylim([0 max(IntM*f)]);        
        set(gca,'xtick',[0,pi/2,pi],...
                'xticklabel',{'0','90','180'});        
        xlabel('$\theta [deg]$','Interpreter','Latex','fontsize',20);
        ylabel('$G(\theta)$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);
        ylim([0 0.5]);

                
        inset2(f1,f2,0.6,[0.25,0.35]);        
        close(f2);          
        
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
 
    function f = OneOverf(t,lambda)
        f = (lambda*(t.^2-(sin(t)).^2).*(pi-t+sin(t).*cos(t))+((pi-t).^2-(sin(t)).^2).*(t-sin(t).*cos(t)))./ ...
                (2*sin(t).*(lambda^2*(t.^2-(sin(t)).^2)+2*lambda*(t.*(pi-t)+(sin(t)).^2)+(pi-t).^2-(sin(t)).^2));
        f(t==0)  = 0;
        f(t==pi) = 0;
    end

    function h = f_interp(y)
        IP = SL.ComputeInterpolationMatrix(SL.CompSpace(y));
        h  = IP.InterPol*(IntM*f)-z;
    end

    function z = GHR(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        z = real(z);
    end

end