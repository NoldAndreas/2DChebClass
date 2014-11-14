function data = SimulationBall(N1,N2,vext)

    disp('** Simulation Ball **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 = 20;   N2 = 20;
        vext  = @VTest;
        R     = 3;
        theta1 = 0;
        theta2 = pi;
        Origin = [0;0];
        N      = [N1;N2];
    end        
    
    
    SG  = Ball(v2struct(R,theta1,theta2,N));
    Int = SG.ComputeIntegrationVector();    
    SG.ComputeInterpolationMatrix((-1:0.02:1)',(0:0.005:0.5)',true,true);
        
    [V,Vdiff,VInt]   = vext(SG.Pts.y1_kv,SG.Pts.y2_kv);   
%    [VP]             = vext(Interp.pts1,Interp.pts2);
           
    %Check Differentiation
    %vplot   = Interp.InterPol*V;        
    %data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
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
    SG.PlotGridLines();    
    SG.PlotGrid();
        
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest(y1,y2)       
            %V           = cos(y1);  VInt = 0;
            V           = sin(y1).^2.*cos(y2).^2;  VInt = 4/3*pi*R^2;
            
            %dVy1        = -y2.*sin(y1);
            %dVy2        = cos(y1);  
            %ddVy1       = -y2.*cos(y1);
            %ddVy2       = zeros(size(y1));
            %dVy1dVy2    = -sin(y1);

            VDiff       = struct();%'dy1',dVy1,'dy2',dVy2,...
                                %'ddy1',ddVy1,'ddy2',ddVy2,...
                                %'dy1dy2',dVy1dVy2);            
                                        
            %sVInt = pi^2*R^2/2;
            %th          = 2*acos(h/R);
            %VInt        = R^2/2*(th-sin(th));            
            
            %C           = sqrt(R^2-h^2);
            %Int        = VInt + sin(C)*(R^2-h^2-C^2+2)-2*cos(C)*C;
             
            
    end

end