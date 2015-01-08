function [data,res] = SimulationBigSegmentAlternative(N1,N2,vext)

    disp('** Simulation SimulationBigSegment **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  20;   N2 = 20; N = 20;        
        N1T = 20; N2T = 20;
        vext  = @VTest;
        R     = 5;
        h     = 2.;
        Dth   = acos(h/R);
        th1   = -pi/2 + Dth;
        th2   = 3/2*pi - Dth;
        Origin = [0;3];
        
        y1Min = -10; y1Max = 10;
        y2Min = -10; y2Max = 10;
    end        
    
    x1Plot = (-1:0.02:1)';
    x2Plot = (-1:0.02:1)';
    
    boxShape           = v2struct(y1Min,y1Max,y2Min,y2Max);
    boxShape.N         = [N1;N2];
    abox               = Box(boxShape);
    [Pts,Diff,Int,Ind] = abox.ComputeAll();
    Interp             = abox.ComputeInterpolationMatrix(x1Plot,x2Plot);

    [V,Vdiff,VInt]     = vext(Pts.y1_kv,Pts.y2_kv); 
    
    
    
	%***************************************************************
    %*** Simulate Arc ***
    %***************************************************************
    ArcShape       = v2struct(R,th1,th2,Origin,N);    
    
    anarc = Arc(ArcShape);
    anarc.ComputeIntegrationVector();
    IP    = abox.SubShapePts(anarc.Pts);
    Int   = anarc.IntSc*IP;
    scatter(anarc.Pts.y1_kv,anarc.Pts.y2_kv,'.')
    disp(['Error in Integration over Outer Path: ',num2str(Int*V-R*(2*pi-2*Dth))]);
    
    %ArcShape.InOut = 'In';
    ArcShape.th1   = th2;       ArcShape.th2   = 2*pi + th1;
    anarc          = Arc(ArcShape);    
    anarc.ComputeIntegrationVector();
    IP             = abox.SubShapePts(anarc.Pts);
    Int            = anarc.IntSc*IP;  
    scatter(anarc.Pts.y1_kv,anarc.Pts.y2_kv,'.')
    disp(['Error in Integration over Inner Path: ',num2str(Int*V-R*2*Dth)]);
        
    %***************************************************************
    %***************************************************************
    
    shape   = v2struct(Origin,R,h);
    shape.N = [N1T,N2T];
    BSA     = BigSegmentAlternative(shape);
    BSA.ComputeIntegrationVector();
    IP      = abox.SubShapePts(BSA.GetCartPts);    
    Int     = BSA.Int*IP;            
    V                = vext(Pts.y1_kv,Pts.y2_kv);   
    [VPS,Vdiff,VInt] = vext(BSA.Pts.y1_kv,BSA.Pts.y2_kv);
           
    %Check Differentiation
    vplot   = IP*V;
    %data    = displayErrorsPos(PtsS,vplot,VPS,V,Vdiff,DiffS,'cart');    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in alternative Segment Integration: ', num2str(data.Int)]);                
                       
    data.N1 = N1; data.N2 = N2;
    
    
    figure;
    BSA.PlotGridLines();    
    BSA.PlotGrid();
    	
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
%     %******************************************
%     %Test Case for Sphere
%     %******************************************
%     disp('Part of a sphere');
%     shape = struct('h',0.0413,'R',0.5,'N',[40 40],'InOut','Out','sphere',true);
%     BS      = BigSegment(shape);        
%     BS.ComputeIntegrationVector();    
%     figure; scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.'); title('BigSegment');
%     
%     %******************************************
%     %Alternative Test Case for Sphere
%     %******************************************
% 	disp('Part of a sphere - Alternative');
%     shape = struct('h',0.0413,'R',0.5,'N',[40 40],'InOut','Out','sphere',true,'Origin',[0 0]);
%     BSA    = BigSegmentAlternative(shape);        
%     BSA.ComputeIntegrationVector();    
%     figure; scatter(BSA.Pts.y1_kv,BSA.Pts.y2_kv,'.'); title('BigSegmentAlternative');
%     
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest(y1,y2)       
            V           = ones(size(y1));%+sin(y1).*y2;  
            dVy1        = cos(y2);
            dVy2        = -y1.*sin(y2);  
            ddVy1       = zeros(size(y1));
            ddVy2       = -y1.*cos(y2);
            dVy1dVy2    = -sin(y2);

            VDiff       = struct('dy1',dVy1,'dy2',dVy2,...
                                'ddy1',ddVy1,'ddy2',ddVy2,...
                                'dy1dy2',dVy1dVy2);            
                                        
            th          = 2*acos(h/R);
            VInt        = (pi*R^2 - R^2/2*(th-sin(th)));
    end

end