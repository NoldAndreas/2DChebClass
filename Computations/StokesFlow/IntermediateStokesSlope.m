function theta = IntermediateStokesSlope(z,lambda)

    global dirData

    fulldir = [dirData filesep 'StokesFlow' filesep];
    
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
    
    if(boolPlot)
        
        plotOpts = struct('plain',true,'linecolor','b','linewidth',2,'linestyle','-');
        f1 =  figure('color','white','Position',[0 0 800 800]);
        plot(y,GHR_lambdaEta0(y),'b','linewidth',2); hold on;
        plotOpts.linecolor = 'r';        
        SL.plot(GHR_lambdaEta(y,1),plotOpts); hold on; 
        
        plotOpts.linecolor = 'm';
        SL.plot(GHR_lambdaEta(y,10),plotOpts); hold on;         
        plotOpts.linecolor = 'k';
        plotOpts.linestyle = ':';
        SL.plot(y.^3/9,plotOpts);
        ylim([0 max(GHR_lambdaEta0(y))]);
        set(gca,'xtick',[0,pi/2,pi],...
                'xticklabel',{'0','90','180'});        
        xlabel('$\theta [deg]$','Interpreter','Latex','fontsize',20);
        ylabel('$G(\theta)$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);                
        
        plotOpts = struct('plain',true,'linecolor','b','linewidth',2,'linestyle','-');        
        f2 =  figure('color','white','Position',[0 0 800 800]);
        plot(y,GHR_lambdaEta0(y),'b','linewidth',2); hold on;
        plotOpts.linecolor = 'r';
        SL.plot(GHR_lambdaEta(y,1),plotOpts); hold on; 
        plotOpts.linecolor = 'm';
        SL.plot(GHR_lambdaEta(y,10),plotOpts); hold on; 
        plotOpts.linecolor = 'k';
        plotOpts.linestyle = ':';
        SL.plot(y.^3/9,plotOpts); 

        ylim([0 max(GHR_lambdaEta0(y))]);        
        set(gca,'xtick',[0,pi/2,pi],...
                'xticklabel',{'0','90','180'});        
        xlabel('$\theta [deg]$','Interpreter','Latex','fontsize',20);
        ylabel('$G(\theta)$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);
        ylim([0 0.5]);

                
        inset2(f1,f2,0.6,[0.25,0.35]);        
        close(f2);          
        
        print2eps([fulldir,'theta_G'],gcf);  
        saveas(gcf,[fulldir,'theta_G.fig']);    
        
        figure('color','white','Position',[0 0 800 800]);
        plot(exp(GHR_lambdaEta0(y)/Ca),y*180/pi,'b','linewidth',2); hold on;
        plot(exp(y.^3/9/Ca),y*180/pi,'b:','linewidth',2);
        xlim([1 max(exp(GHR_lambdaEta0(y)/Ca))])        
        xlabel('$r/\hat L$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta(r)[^\circ]$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',15);
        
        print2eps([fulldir,'r_theta'],gcf);
        saveas(gcf,[fulldir,'r_theta.fig']);
        %subplot(1,2,2); plot(theta,z,'ko','markersize',6,'MarkerFaceColor','r');
    end
end