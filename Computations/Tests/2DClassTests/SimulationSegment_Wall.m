function data = SimulationSegment_Wall(N1,N2,vext)

    disp('** Simulation Segment M1 **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  20;   N2 = 20;
        vext  = @VTest;
        R     = 3;
        h     = -1; 
        Origin = [0;0];
        N      = [N1;N2];
    end            
    
    SG                 = Segment(v2struct(Origin,R,h,N));
    [Pts,Diff,Int,Ind] = SG.ComputeAll();
    Interp             = SG.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
        
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv);
    [VP]             = vext(Interp.pts1,Interp.pts2);
           
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
        
    subplot(2,2,1); SG.doPlots(Vdiff.dy1,'SC');
    subplot(2,2,2); SG.doPlots(Diff.Dy1*V,'SC');
    subplot(2,2,3); SG.doPlots(Diff.Dy1*V-Vdiff.dy1,'SC');    
    subplot(2,2,4); SG.doPlots(V,'SC');
               
    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color            
    SG.doPlots(V,'SC');    
    title('Interpolation');      
    figure
    SG.PlotGridLines(); %scatter(SG.Pts.y1_kv,SG.Pts.y2_kv,'.')
    
    %***************************************************************
    %   Check integration of spherical cap:
    %***************************************************************             
    shape  = v2struct(Origin,R,h,N);
    shape.sphere = true;
    SG_Sph = Segment(shape);
    intSph = SG_Sph.ComputeIntegrationVector();    
    SG_Sph.PlotGrid();
    
    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest(y1,y2)       
            V           = ones(size(y1))+y2.*cos(y1);  
            dVy1        = -y2.*sin(y1);
            dVy2        = cos(y1);  
            ddVy1       = -y2.*cos(y1);
            ddVy2       = zeros(size(y1));
            dVy1dVy2    = -sin(y1);

            VDiff       = struct('dy1',dVy1,'dy2',dVy2,...
                                'ddy1',ddVy1,'ddy2',ddVy2,...
                                'dy1dy2',dVy1dVy2);            
                            
            disp('Integration only works if baseline = h');
            th          = 2*acos(h/R);
            VInt        = R^2/2*(th-sin(th));            
            
            C           = sqrt(R^2-h^2);
            VInt        = VInt + sin(C)*(R^2-h^2-C^2+2)-2*cos(C)*C;
             
            
    end

end