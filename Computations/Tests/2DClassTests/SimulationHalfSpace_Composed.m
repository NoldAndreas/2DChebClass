function data = SimulationHalfSpace_Composed()%N1,N2,L1,L2,vext)

    disp('** Simulation Half Space Composed**');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)
        N      = [40;30];
        N2bound = 10;
        L1 = 2;    L2 = 1;        
        y2Min = 0; h = 1;        
        vext  = @Vext15;       
        alpha = pi/4;
    end                
    
    HSC         = HalfSpace_Composed(v2struct(L1,L2,N,N2bound,y2Min,h,alpha));            
    HSCCartPts  = HSC.GetCartPts();
    V           = vext(HSCCartPts.y1_kv,HSCCartPts.y2_kv);
        
    HSC.plot(V,'SC');%,'SC'); 
    
    PhysAreaBX  = struct('y1Min',-3*L1,'y1Max',3*L1,...
                        'y2Min',0,'y2Max',10,...
                        'N',[40,40]);
    BX          = Box(PhysAreaBX);
    BXCartPts   = BX.GetCartPts();    
    
    
    Interp      = HSC.SubShapePtsCart(BX.Pts);    
    VP          = vext(BXCartPts.y1_kv,BXCartPts.y2_kv);

    vplot       = Interp*V;    
    data        = displayErrorsPos(HSC.Pts,vplot,VP,V);
    
end