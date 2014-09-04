function res = LubricationForcedWettingInnerFinite(delta)


    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'ForcedWetting'],'ORG');    
    close all;
    
    %% Parameters
    if(nargin == 0)
        delta   = 0.1;
    end
    N       = 500;
    
    IntM    = 0;    
    Dy = 0; DDy = 0;
    
    y= 0;
    h1 = 0; h2 = 0;
    
    SolveInner_ThirdOrder_L(100);    
    subplot(1,2,1);
    SolveInnerL(100);
    
    [res.L20.y,res.L20.hP] = SolveInnerL(20);
    [res.L50.y,res.L50.hP] = SolveInnerL(50);
    [res.L100.y,res.L100.hP] = SolveInnerL(100);
    [res.L500.y,res.L500.hP] = SolveInnerL(500);
    [res.L5000.y,res.L5000.hP] = SolveInnerL(5000);    
    
     function SolveInner_ThirdOrder_L(L)

        shape      = struct('N',N,'yMin',0,'yMax',L);    
        plotShape  = struct('yMin',0,'yMax',L,'N',N);

        SL             = SpectralLine(shape);    
        [Pts,Diff,Int] = SL.ComputeAll(plotShape);                
        
        IntM = zeros(length(y));
        vh     = zeros(1,length(y));
        for i = 2:length(y)
            %Matrix for integral int(v(y),y=y..Y)
            hh          = y(i) - y(i-1);
            vh([i-1,i]) = vh([i-1,i]) + hh/2;
            IntM(i,:)    = vh;
        end    

        y              = Pts.y; 
        Dy             = Diff.Dy;
        DDy            = Diff.DDy;              
        
        h0 = y;
        h1 = -(1/2)*y.^2.*log(y)+(1/2)*log(y+1).*(y+1).^2-(1/2)*y;
        h1(1) = 0;
        h2 = -(1/2*(y+1)).*((y+2).*log(y)+2*y-1).*log(y+1)+(1/6*(-6*y.^2-6*y-6)).*dilog(y+1)-(1/6)*y.*((-6*y-6).*log(y)+y*pi^2+9);
        h2(1) = 0; 
        
        
        plot(y,h1.*(2*y+1)./((y.*(y+1)).^2)); hold on;
        plot(y,((-3*y.^2-3*y-1).*h1.^2+(2*(y+1)).*y.*(y+1/2).*h2)./(y.^3.*(y+1).^3)); 
        
        h2P  = fsolve(@ODE2,zeros(N,1));
        
       % h2 = IntM*h2P;
        
        h3P  = fsolve(@ODE3,zeros(N,1));
        
        SL.doPlots(h2P);
        SL.doPlots(Dy*h2);
        
        subplot(1,2,1);
        plot(y,Dy*h0,'-'); hold on;
        plot(y,Dy*(h0+delta*h1),'--');
        plot(y,Dy*(h0+delta*h1+delta^2*h2),':');
        plot(y,Dy*(h0+delta*h1)+delta^2*h2P,':m');
        plot(y,Dy*(h0+delta*h1+delta^2*h2)+delta^3*h3P,'-.');                       
        
        yP = (1:0.01:L)';
        plot(yP,(3*delta*log(yP)+1+3*delta).^(1/3),'r');
        subplot(1,2,2);
        hold on; SL.doPlots(h3P);
    
 
    end
    
    function [y,hP] = SolveInnerL(L)
        shape      = struct('N',N,'yMin',0,'yMax',L);    
        plotShape  = struct('yMin',0,'yMax',L,'N',N);

        SL             = SpectralLine(shape);    
        [Pts,Diff,Int] = SL.ComputeAll(plotShape);                

        y              = Pts.y; 
        Dy             = Diff.Dy;
        DDy            = Diff.DDy;      

        IntM = zeros(length(y));
        vh     = zeros(1,length(y));
        for i = 2:length(y)
            %Matrix for integral int(v(y),y=y..Y)
            hh          = y(i) - y(i-1);
            vh([i-1,i]) = vh([i-1,i]) + hh/2;
            IntM(i,:)    = vh;
        end    

        hP  = fsolve(@ODE,ones(N,1));
        SL.doPlots(hP);
    end

    hold on;
    y = (1:0.01:5000)';
    plot(y,(3*delta*log(y)+1+3*delta).^(1/3));
    
    function z = ODE3(hP)
        z      = ((-3*y.^2-3*y-1).*h1.^2+(2*(y+1)).*y.*(y+1/2).*h2)./(y.^3.*(y+1).^3) - DDy*hP; 
        z(1)   = hP(1);
        z(end) = Dy(end,:)*hP;
    end
    function z = ODE2(hP)
        z      = h1.*(2*y+1)./((y.*(y+1)).^2) - DDy*hP; 
        z(1)   = hP(1);
        z(end) = Dy(end,:)*hP;
    end    
    function y = ODE(hP)
        
        h = IntM*hP;        
        y = delta./(h.^2+h) + DDy*hP;
        
        % Boundary conditions 
        y(1)   = hP(1)-1;
        y(end) = Dy(end,:)*hP;
        
    end

    

end