function data = SimulationBigSegment_WedgeAndTriangle()

    close all;
    N1T = 10; N2T = 10;
    
    
	shape.Origin = [3,-4];
    shape.R      = 2;
    shape.th1 = 60*pi/180;
    shape.th2 = 30*pi/180;
    shape.NW = [4*N1T,N2T];   
    shape.NT = [N1T,N2T];   
    shape.sphere = true;
    
    BS      = BigSegment(shape);    
    figure;    
    scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.'); title('Composed Map');
       
    figure;    
    
    subplot(2,2,1);
    shape.Wall_Y = -5;
    shape.Wall_VertHor = 'horizontal';    
    BS = BigSegment(shape);
    scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.'); 
    
    subplot(2,2,2);
    shape.Wall_Y = -2.5;
    shape.Wall_VertHor = 'horizontal';    
    BS = BigSegment(shape);
    scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.'); 
    
    subplot(2,2,3);
    shape.Wall_Y = 2;
    shape.Wall_VertHor = 'vertical';    
    BS = BigSegment(shape);
    scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.'); 
    
    subplot(2,2,4);
    shape.Wall_Y = 4.5;
    shape.Wall_VertHor = 'vertical';    
    BS = BigSegment(shape);
    scatter(BS.Pts.y1_kv,BS.Pts.y2_kv,'.');     
    
    
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