

    %Initialize half space
    Phys_Area = struct('L1',3,'L2',2,'N',[20;20]);

    PlotArea       = struct('y1Min',-5,'y1Max',5,...
                            'N2',100,'N1',100,'y2Min',0,'y2Max',4);        
    HS                             = HalfSpace(Phys_Area);
   % [Pts,Diff,Int,Ind,Interp,Conv] = HS.ComputeAll(PlotArea,@f2);
    
    %Initialize infinite annulus
	Origin       = [0,3];
    N            = [15,16];
    sphere       = true;
    shapeDC      = v2struct(Origin,N,sphere,L);
    shapeDC.RMin = 1;
    DC           = InfAnnulus(shapeDC);   

    hold on
    area = Intersect(HS,DC);
    
    area.shape.PlotGridLines();
    xlim([-5 5]); ylim([0 5]);
    
    
    gridy1 = ((Origin(1)-5):0.1:(Origin(1)+5))';
    gridy2 = (0:0.1:(Origin(2)+5))';
    Grid.y1_kv = kronecker(gridy1,ones(size(gridy2)));
    Grid.y2_kv = kronecker(ones(size(gridy1)),gridy2);
    
    
    
        