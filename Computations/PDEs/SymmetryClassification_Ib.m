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
    lambda1 = 3;
    lambda2 = 1;
    
    
    %% Initialization
    %     
        
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
end

