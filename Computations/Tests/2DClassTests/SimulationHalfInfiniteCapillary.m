function SimulationHalfInfiniteCapillary

    disp('** Simulation Half Infinite Capillary **');

    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
            
    vext  = @Vext7;
    
	Phys_Area = struct('L1',1,'N',[40;30],...
                       'y2Min',0,'y2Max',4);
    PlotArea  = struct('y1Min',0,'y1Max',5,...
                       'N2',100,'N1',100,'y2Min',0,'y2Max',4);    
    
    HIC                       = HalfInfCapillary(Phys_Area);
    [Pts,Diff,Int,Ind,Interp] = HIC.ComputeAll(PlotArea);        
    
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv);%,intBound,'cart');    
    [VP]             = vext(Interp.pts1,Interp.pts2);           
                   
    %Check Differentiation
    vplot     = Interp.InterPol*V;        
    displayErrors(vplot,VP,V,Vdiff,Diff,'cart');
    
    %Check Interpolation    
    HIC.doPlots(V,'SC');                    
    
    %Check Integration    
    display([' Error in Integration: ', num2str(Int(Int < inf)*V(Int < inf)-VInt)]);                
            
                                                                     
 
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f2(x,y)
        z = exp(-(x.^2+y.^2));
    end
end