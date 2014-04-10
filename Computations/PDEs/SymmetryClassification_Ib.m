function SymmetryClassification_Ib()
%%  Solution I.b for a linear shear flow

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
%

    global dirData
    AddPaths();
    ChangeDirData([dirData filesep 'SymmetryClassification'],'ORG');

    %% Parameters
    lambda1 = 0;    lambda2 = 1;  omega   = -1i;    
    
    C1 = 1;     C2 = 0;
        
    x0 = 0;     y0 = 0; 
    AL = 1;     VL = 0; %Parameters for base flow U(y) = AL*y + VL
    
    %% Initialization
    %     
        
    shapeBox = struct('y1Min',-5,'y1Max',5,'N',[80,50],...
                      'y2Min',0.1,'y2Max',1);
    plotBox  = shapeBox;
    plotBox.N = [80,80];
    BX        = Box(shapeBox);
    Pts       = BX.ComputeAll(plotBox);        
    
    t = [0,0.5,1,1.5,2];
    
    for i = 1:length(t)
        subplot(2,2,i);
        p = Psi(Pts.y1_kv,Pts.y2_kv,t(i));
        BX.doPlots(real(p),'contour',struct('clabel',true));    pbaspect([1 1 1]);  
        %subplot(1,2,2);     BX.doPlots(imag(p),'contour',struct('clabel',true));    pbaspect([1 1 1]);
    end
    
    %% ****************************************************************
    shapeBox = struct('y1Min',-5,'y1Max',5,'N',[50,70],...
                      'y2Min',0.0,'y2Max',10);
    plotBox  = shapeBox;
    plotBox.N = [80,80];
                                    
    %BX       = Box(shapeBox);    
    %[PtsBx]  = BX.ComputeAll(plotBox);
    
    shapeHS  = struct('L1',3,'L2',3,'N',[50,70]);
    HS       = HalfSpace(shapeHS);
    [PtsHS]  = HS.ComputeAll(plotBox);
    
    y1       = PtsHS.y1_kv;
    y2       = PtsHS.y2_kv;
    
    F = RHS(y1,y2);
    HS.doPlots(F,'SC');
    
   
          
    %% Right hand side of ODE
    function h = RHS(x,y)
        z = x./y + log(y)/lambda2;                 
        
        z((y==0) & (x>0))  = inf;        
        z((y==0) & (x<=0)) = -inf;        
        z(y == inf)        = -inf/lambda2;
        
        h = (z.^2).*exp(-z.^2).*(y.^(lambda1/lambda2-2));
        h(y==inf) = 0;        
        h(z==inf) = 0;        
        h(z==-inf) = 0;
    %    h(y==0) = 0;
    end

    function p = Psi(x,y,t)
        z1   = lambda2*(x-x0-(y0*AL+VL)*t)./y;
        eta = log(y-y0)-lambda2*AL*t;        
        
        z2  = (-1i*omega)^(1/3)*(1i/(4*omega) + 1i*omega*lambda2^2 - 1 - eta);
        
        p   = C1*exp(1i*lambda2*(z1+lambda1*AL*t)+eta/2.*(1-1i*eta*omega)).*airy(0,z2) + ...
              C2*exp(1i*lambda2*(z1+lambda1*AL*t)+eta/2.*(1-1i*eta*omega)).*airy(2,z2);
    end
end

