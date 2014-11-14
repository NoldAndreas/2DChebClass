function data = SimulationSegment2(N1,N2,vext)

    disp('** Simulation Sehment 2 **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  30;   N2 = 30;
        vext  = @VTest;
        R     = 3;
        theta1 = 0;
        theta2 = pi*3/4;
        Origin = [0;0];
        N      = [N1;N2];
    end        
    
    
    SG  = Segment2(v2struct(R,theta1,theta2,N));
    Int = SG.ComputeIntegrationVector();
    SG.ComputeInterpolationMatrix((-1:0.02:1)',(0:0.01:0.5)',true,true);
    
    [V,Vdiff,VInt]   = vext(SG.Pts.y1_kv,SG.Pts.y2_kv);
    PtsCart          = SG.GetCartPts();
    VCart            = vext(PtsCart.y1_kv,PtsCart.y2_kv);
%    [VP]             = vext(Interp.pts1,Interp.pts2);
           
    %Check Differentiation
    %vplot   = Interp.InterPol*V;        
    %data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Integration
    data.Int = abs(Int*VCart-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                            
    
%    subplot(2,2,1); doPlots_SC(Interp,Pts,Vdiff.dy1);
%    subplot(2,2,2); doPlots_SC(Interp,Pts,Diff.Dy1*V);
%    subplot(2,2,3); doPlots_SC(Interp,Pts,Diff.Dy1*V-Vdiff.dy1);    
               
    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color            
    SG.plot(V,true);        
    title('Interpolation');      
    
    figure 
    SG.PlotGridLines(); hold on;
    SG.PlotGrid();
    
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest(y1,y2)       
            V           = cos(y2).*y1;  
            dVy1        = -y2.*sin(y1);
            dVy2        = cos(y1);  
            ddVy1       = -y2.*cos(y1);
            ddVy2       = zeros(size(y1));
            dVy1dVy2    = -sin(y1);

            VDiff       = struct('dy1',dVy1,'dy2',dVy2,...
                                'ddy1',ddVy1,'ddy2',ddVy2,...
                                'dy1dy2',dVy1dVy2);            
                            
            disp('Integration only works if baseline = h');
            VInt = 0;
            %th          = 2*acos(h/R);
            %VInt        = R^2/2*(th-sin(th));            
            
            %C           = sqrt(R^2-h^2);
            %Int        = VInt + sin(C)*(R^2-h^2-C^2+2)-2*cos(C)*C;
             
            
    end

end