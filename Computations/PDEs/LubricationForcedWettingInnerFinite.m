function LubricationForcedWettingInnerFinite


    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'ForcedWetting'],'ORG');    
    close all;
    
    %% Parameters
    delta   = 0.3;
    N       = 300;
    
    IntM    = 0;    
    Dy = 0; DDy = 0;
    
    SolveInner_ThirdOrder_L(20);
    
    SolveInnerL(10);
    SolveInnerL(20);
    SolveInnerL(50);
    SolveInnerL(100);
    SolveInnerL(500);
    SolveInnerL(5000);    
    
    function [y,hP] = SolveInner_ThirdOrder_L(L)

        shape      = struct('N',N,'yMin',0,'yMax',L);    
        plotShape  = struct('yMin',0,'yMax',L,'N',N);

        SL             = SpectralLine(shape);    
        [Pts,Diff,Int] = SL.ComputeAll(plotShape);                

        y              = Pts.y; 
        Dy             = Diff.Dy;
        DDy            = Diff.DDy;              
        
        h0 = y;
        h1 = -(1/2)*y.^2.*log(y)+(1/2)*log(y+1).*(y+1).^2-(1/2)*y;
        h1(1) = 0;
        h2 = -(1/2*(y+1)).*((y+2).*log(y)+2*y-1).*log(y+1)+(1/6*(-6*y.^2-6*y-6)).*dilog(y+1)-(1/6)*y.*((-6*y-6).*log(y)+y*pi^2+9);
        h2(1) = 0;
        Dh3y3 = ((-3*y.^2-3*y-1).*h1.^2+(2*(y+1)).*y.*(y+1/2).*h2)./(y.^3.*(y+1).^3); 
        
        plot(y,Dy*h0,'-'); hold on;
        plot(y,Dy*(h0+h1),'--');
        plot(y,Dy*(h0+h1+h2),':');
        plot(y(2:end),Dy(2:end,:)*(h0+h1+h2)+Dh3y3(2:end),'-.');
        
        yP = (1:0.01:L)';
        plot(yP,(3*delta*log(yP)+1+3*delta).^(1/3),'r');
    
 
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
    
    function y = ODE(hP)
        
        h = IntM*hP;        
        y = delta./(h.^2+h) + DDy*hP;
        
        % Boundary conditions 
        y(1)   = hP(1)-1;
        y(end) = Dy(end,:)*hP;
        
    end

    

end