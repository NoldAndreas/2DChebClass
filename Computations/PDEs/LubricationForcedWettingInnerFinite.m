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
    
    SolveInnerL(10);
    SolveInnerL(20);
    SolveInnerL(50);
    SolveInnerL(100);
    SolveInnerL(500);
    
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
    
    function y = ODE(hP)
        
        h = IntM*hP;        
        y = delta./(h.^2+h) + DDy*hP;
        
        % Boundary conditions 
        y(1)   = hP(1)-1;
        y(end) = Dy(end,:)*hP;
        
    end

end