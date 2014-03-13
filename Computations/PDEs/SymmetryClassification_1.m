function SymmetryClassification_1()
%Solution I.a for a linear shear flow

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
% 


    %% Parameters
    k = 0.5;
    N = 100;
    L = 4;       
    

    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape     = struct('N',N,'L',L);    
    plotShape = struct('yMin',-5,'yMax',5,'N',200);

    IS         = InfSpectralLine(shape);    
    [Pts,Diff] = IS.ComputeAll(plotShape);        
    
    y    = Pts.y(2:end-1);    
    Dy   = Diff.Dy(2:end-1,2:end-1);
    DDy  = Diff.DDy(2:end-1,2:end-1);        
        
    shapeBox = struct('y1Min',-5,'y1Max',5,'N',[40,40],...
                      'y2Min',0,'y2Max',5);
    plotBox = struct('y1Min',-5,'y1Max',5,'N',[100,100],...
                     'y2Min',0,'y2Max',5);                  
                  
    BX       = Box(shapeBox);
    [PtsBx]  = BX.ComputeAll(plotBox);
    
    zeta     = PtsBx.y1_kv./PtsBx.y2_kv;
    IP       = IS.InterpolationMatrix_Pointwise(zeta);
    
    %% Solve ODE
    %
    % # set up operator 
    % # in for-loop, set up and right hand side for time t
    % # compute and plot solution
    
    A = k*(k-1)*diag(1./(1+y.^2))+2*(1-k)*diag(y./(1+y.^2))*Dy+DDy;
    
    for t = 0:0.1:3
        h = diag(1./(1+y.^2))*RHS(Pts.y(2:end-1)-t);
        f = A\h;    
        f = [0;f;0];
        
        %compute 
        Psi = (PtsBx.y2_kv.^k).*(IP*f);
        
        subplot(1,2,1);
        hold off;       
        IS.doPlots(f);
        
        subplot(1,2,2);
        BX.doPlots(Psi,'contour',struct('clabel',true));
        
        pause(0.1);
    end
    
      
    

    %% Right hand side of ODE
    function h = RHS(y)
        h = exp(-2*y.^2);
    end
end

