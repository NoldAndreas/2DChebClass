function SolveLubrication()

    N    = 60;
    L    = 4;
    hMax = 6;
    a    = 2;
    
    xxp = (-1:0.01:1)';    
    xp  = ClenCurtFlip(N);
    IP  = barychebevalMatrix(xp,xxp);
    yp  = xp*L;
    yyp = xxp*L;
    
    yp2  = (1+xp)*a*L + L;
    yyp2 = (1+xxp)*a*L + L;
    
    Diff = barychebdiff(xp);
    Dy   = Diff.Dx/L;
    Dy2  = Diff.Dx/(a*L);
    
    h0      = fsolve(@DisjoiningPressure,0.9);
    disp(['h0 = ',num2str(h0)]);
    theta_0 = acos(1+BindingPotential(h0));
    disp(['Equilibrium Contact Angle: ',num2str(theta_0*180/pi),' [deg]']);
    
    opts = optimset('MaxFunEvals',200000,'MaxIter',1000);
    h1 = fsolve(@f,log(hMax*ones(length(yp)-1,1)-h0));    
    h1 = [h0+exp(h1);hMax];
    
    h2 = fsolve(@f2,log(hMax*ones(length(yp)-1,1)-h0));
    h2 = [hMax;h0+exp(h2)];
    
    %h = fsolve(@f,hMax*ones(length(yp)-1,1),opts);    
    h = [h1;h2];
    
    figure('color','white','Position',[0 0 600 600]);
    hIDL = (0.8:0.01:20);
    plot(hIDL,BindingPotential(hIDL),'linewidth',2); hold on;
    plot([0 max(hIDL)],[0 0],'k--','linewidth',1.5);
    set(gca,'fontsize',20);
    ylim([-0.3 0.35]);
    %set(gca,'YTick',-0.1:0.1:0.3);
    xlabel('${h}/h_0$','Interpreter','Latex','fontsize',25);
    ylabel('${w}/{\gamma_{{lv}}}$','Interpreter','Latex','fontsize',35);
    xlim([0 20]);
    
	print2eps('DisjoiningPotential',gcf);
	saveas(gcf,'DisjoiningPotential.fig');
    
    figure('color','white','Position',[0 0 600 600]);
%    plot([yp;yp2],h,'o'); hold on;
    yG = [yyp;yyp2];
    plot(yG,[IP*h1;IP*h2]/h0,'linewidth',2); hold on;
    %plot(yyp,IP*h);
    
    plot(yG,(max(h)+tan(theta_0)*(yG-max(yG)))/h0,'k--','linewidth',2);
    plot([min(yyp) max(yyp2)],[h0 h0]/h0,'m--','linewidth',1.5)
    set(gca,'fontsize',20);
    xlim([min(yyp) max(yyp2)]);
    xlabel('$x$','Interpreter','Latex','fontsize',25);
    ylabel('$h/h_0$','Interpreter','Latex','fontsize',25);
        
	print2eps('FilmThickness',gcf);
	saveas(gcf,'FilmThickness.fig');

    function z = f(h)        
        h = [h0+exp(h);hMax];
        z = Dy*h - dh(h);
        z = z(1:end-1);
    end

    function z = f2(h)
        h = [hMax;h0+exp(h)];
        z = Dy2*h - dh(h);
        z = z(2:end);
    end

    
    function z = dh(h)        
        z = sqrt( (cos(theta_0) - BindingPotential(h)).^(-2) - 1 );
    end
    
    
    function [w,dw] = BindingPotential(h)        
        
        w   = 7./(18*h.^8) - 2./h.^2 + 12./(h+2).^2;
        dw  = -8*7./(18*h.^9) + 4./h.^3 - 24./(h+2).^3;
        %a  = 1.05;
        %w  = 1./h.^8 - 2*a./h.^2 + 12./(h+2).^2;
        %dw = -8./h.^9 + 4*a./h.^3 - 24./(h+2).^3;
    end

    function p = DisjoiningPressure(h)        
        [h1s,p] = BindingPotential(h);
    end
    
end

