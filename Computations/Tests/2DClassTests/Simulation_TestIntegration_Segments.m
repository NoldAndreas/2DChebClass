function Simulation_TestIntegration_Segments

    Wall_VertHor = 'horizontal';
    Wall_Y = -5;
    R      = 2;

    if(strcmp(Wall_VertHor,'horizontal'))
        plot([-5 5],[Wall_Y,Wall_Y]);   hold on; 
        plot([-5 5],[Wall_Y+R,Wall_Y+R],'r--','LineWidth',2);    
        plot([-5 5],[Wall_Y-R,Wall_Y-R],'r--','LineWidth',2);    
    
        xlim([-3,3])
        ylim([Wall_Y-2*R,Wall_Y+2*R])
        axis equal
    elseif(strcmp(Wall_VertHor,'vertical'))
        plot([Wall_Y,Wall_Y],[-5 5]);   hold on; 
        plot([Wall_Y+R,Wall_Y+R],[-5 5],'r--','LineWidth',2);    
        plot([Wall_Y-R,Wall_Y-R],[-5 5],'r--','LineWidth',2);    
    
        ylim([-3,3])
        xlim([Wall_Y-2*R,Wall_Y+2*R])
        axis equal        
    end
    
    
    
    while(true)
        
        [y10,y20] = ginput(1);
    
        shape = struct('Origin',[y10,y20],'R',R,'N',[20,20],...
                        'Wall_VertHor',Wall_VertHor,'Wall_Y',Wall_Y);

        SW  = SegmentWall(shape);
        SW.ComputeIntegrationVector();
        
        shape.NW = [20,20];
        shape.NT = [20,20];
        BSW = BigSegment(shape);
        BSW.ComputeIntegrationVector();

        SW.PlotGrid(); hold on;
        BSW.PlotGrid();
        
        %Compute integration 
        SW_CartPts  = SW.GetCartPts();
        BSW_CartPts = BSW.GetCartPts();
        
        V_SW  = VTest1(SW_CartPts.y1_kv - y10,SW_CartPts.y2_kv - y20);
        V_BSW = VTest1(BSW_CartPts.y1_kv - y10,BSW_CartPts.y2_kv - y20);
        [~,VInt]  = VTest1([],[]);
        disp(['Error for VTest1: ',num2str(SW.Int*V_SW + BSW.Int*V_BSW - VInt)]);
        
        V_SW  = VTest2(SW_CartPts.y1_kv - y10,SW_CartPts.y2_kv - y20);
        V_BSW = VTest2(BSW_CartPts.y1_kv - y10,BSW_CartPts.y2_kv - y20);
        [~,VInt]  = VTest2([],[]);
        disp(['Error for VTest2: ',num2str(SW.Int*V_SW + BSW.Int*V_BSW - VInt)]);
        
        V_SW  = VTest3(SW_CartPts.y1_kv - y10,SW_CartPts.y2_kv - y20);
        V_BSW = VTest3(BSW_CartPts.y1_kv - y10,BSW_CartPts.y2_kv - y20);
        [~,VInt]  = VTest3([],[]);
        disp(['Error for VTest3: ',num2str(SW.Int*V_SW + BSW.Int*V_BSW - VInt)]);
        
        V_SW  = VTest4(SW_CartPts.y1_kv - y10,SW_CartPts.y2_kv - y20);
        V_BSW = VTest4(BSW_CartPts.y1_kv - y10,BSW_CartPts.y2_kv - y20);
        [~,VInt]  = VTest4([],[]);
        disp(['Error for VTest4: ',num2str(SW.Int*V_SW + BSW.Int*V_BSW - VInt)]);
    end
    
    
    
    function [V,VInt] = VTest1(y1,y2)
        V       = y1.^2 + y2.^2;
        VInt    = pi*R^4/2;
    end
    
	function [V,VInt] = VTest2(y1,y2)
        V       = y1.^2;
        VInt    = pi*R^4/4;
    end

    function [V,VInt] = VTest3(y1,y2)
        V       = y2.*sin(y1);
        VInt    = 0;
    end

    function [V,VInt] = VTest4(y1,y2)
        V       = y1;
        VInt    = 0;
    end

end