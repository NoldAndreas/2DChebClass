function SymmetryClassification_Ib_Eigenfunctions()
%%  Eigensolution I.b for a linear shear flow

%% ODE to solve
% 
% $$(x^2+1)\frac{\partial^2}{\partial x^2} + \frac{\partial}{\partial y} + \frac{\partial^2}{\partial y^2} + 2x \left(\frac{\partial^2}{\partial x\partial y} + \frac{\partial}{\partial x}\right) = h(y-\lambda_2x)e^{-cx}$$
%

    global dirData
    AddPaths();
    ChangeDirData([dirData filesep 'SymmetryClassification'],'ORG');

    %% Parameters
    c = 0;    lambda2 = 1;  
    
    %% Initialization
    %     
    PhysArea = struct('N',[25;25],'L1',2,'L2',2.,'y2wall',0.,...
                       'N2bound',16,'h',1,'L2_AD',1.);

    PlotArea = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);

	IS                 = InfSpace(PhysArea);
    [Pts,Diff,Int,Ind] = IS.ComputeAll(PlotArea);    
    
    IS.doPlots(rhs(Pts.y1_kv,Pts.y2_kv),'contour');
    
    
          
    %% Right hand side of ODE
    function z = rhs(y1,y2)
        z            = h(y2-lambda2*y1).*exp(-c*y1);        
        z(Ind.bound) = 0;
    end

    function z = h(dy)
        z = exp(-(dy/3).^2);
    end
end

