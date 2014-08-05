function theta =IntermediateStokesSlope(z,lambda,thetam)
    if(nargin == 0)
        z      = (0.01:0.02:0.2)';
        lambda = 1;
        boolPlot = true;
    else
        boolPlot = false;
    end
        
    N          = 100;
    shape      = struct('N',N,'yMin',0,'yMax',pi);    
    plotShape  = struct('yMin',0,'yMax',pi,'N',200);

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
    
    IP_m  = SL.ComputeInterpolationMatrix(SL.CompSpace(thetam));
    z     = z + IP_m.InterPol*f;
    theta = fsolve(@f_interp,(9*z).^(1/3));
    
    if(boolPlot)
        subplot(1,2,1); SL.doPlots(f); hold on; SL.doPlots(y.^2/3);
        subplot(1,2,2); SL.doPlots(IntM*f);  hold on; SL.doPlots(y.^3/9); 

        ylim([0 max(IntM*f)]);
        xlabel('$\theta$','Interpreter','Latex');
        ylabel('$G(\theta)$','Interpreter','Latex');

        subplot(1,2,2); plot(theta,z,'ko','markersize',6,'MarkerFaceColor','r');
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