function data = SimulationHalfSpace_Composed(N1,N2,L1,L2,vext)

    disp('** Simulation Half Space Composed**');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)
        N      = [30;20];
        N2bound = 10;
        L1 = 2;    L2 = 1;        
        y2Min = 0; h = 1;        
        vext  = @Vext15;        
    end                
    
    HSC = HalfSpace_Composed(v2struct(L1,L2,N,N2bound,y2Min,h));        
    V   = vext(HSC.Pts.y1_kv,HSC.Pts.y2_kv); 
    [VIP,ptsPlot] = HSC.doPlots(V,'SC');        
        
    VP  = vext(ptsPlot.y1_kv,ptsPlot.y2_kv); 
    
    
    [data.InterPol,iMaxErr] = max(abs(VIP - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol) ,...
             ' at y1 = ',num2str(ptsPlot.y1_kv(iMaxErr)),...
             ' at y2 = ',num2str(ptsPlot.y2_kv(iMaxErr))]);        
    
end