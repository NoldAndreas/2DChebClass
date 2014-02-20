function data = SimulationHalfSpaceSkewed_Composed()

    disp('** Simulation Half Space Skewed Composed**');
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
        alpha = pi/4;
    end                
    
    HSC = HalfSpaceSkewed_Composed(v2struct(L1,L2,N,N2bound,y2Min,h,alpha));        
    V   = vext(HSC.Pts.y1_kv,HSC.Pts.y2_kv); 
    HSC.doPlots(V);%,'SC');        
    
    
    y1_I = (-5:0.1:5)';
    y2_I = (0:0.1:5)';
    ptsI.y1_kv  = kron(y1_I,ones(size(y2_I)));
    ptsI.y2_kv  = kron(ones(size(y1_I)),y2_I);    
        
    VP  = vext(ptsI.y1_kv,ptsI.y2_kv); 
    IP  = HSC.SubShapePts(ptsI);
    
    
    [data.InterPol,iMaxErr] = max(abs(IP*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol) ,...
             ' at y1 = ',num2str(ptsI.y1_kv(iMaxErr)),...
             ' at y2 = ',num2str(ptsI.y2_kv(iMaxErr))]);        
    
end