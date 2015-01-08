function SymmetryClassification_Ib_Homogeneous()
%%  Solution I.b for a linear shear flow

%% ODE to solve
% 
% $$k(k-1)f+2y(1-k)f'+(1-y^2)f''=h(y-t)$$
%

    global dirData
    AddPaths();
    ChangeDirData([dirData filesep 'SymmetryClassification'],'ORG');

    %% Parameters
    k = 0;
    %% Initialization
    %     
        
    shapeBox = struct('y1Min',-0.5,'y1Max',0.5,'N',[80,50],...
                      'y2Min',0.,'y2Max',1);
    plotBox  = shapeBox;
    plotBox.N = [100,100];
    BX        = Box(shapeBox);
    Pts       = BX.ComputeAll(plotBox);        
    
    figure('Position',[0 0 800 800],'color','white');
    p = Psi(Pts.y1_kv,Pts.y2_kv);
    BX.plot(p,'contour',struct('clabel',true,'linecolor','k'));    pbaspect([1 1 1]);  
   
    xlabel('$x$','Interpreter','Latex','fontsize',25);
    ylabel('$y$','Interpreter','Latex','fontsize',25);   
    
    print2eps([dirData filesep 'SelfSimilarSolution_Ib_Hom'],f1);
    saveas(f1,[dirData filesep 'SelfSimilarSolution_Ib_Hom.fig']);
    
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
    HS.plot(F,'SC');
    
   
          
    %% Right hand side of ODE
    function p = Psi(x,y)        
        p = atan(x./y);        
    end
end

