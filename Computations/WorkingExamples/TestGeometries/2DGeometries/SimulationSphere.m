function [data,res] = SimulationSphere()

    disp('** Simulation Ball - Sphere**');
    AddPaths();
    
    close all;
    
    %Initialization
    N1 = 20;
    N2 = 20;
    vext  = @VTestt;
    R      = 1;
    theta1 = pi/5;%pi/2;
    theta2 = pi*3/4;
    Origin = [0;0];
    N      = [N1;N2];    
    volume = false;
    
    shape  = v2struct(R,theta1,theta2,N,volume);
    
    SG  = Sphere(shape);
    Pts = SG.Pts;
    Int = SG.ComputeIntegrationVector();    
    Interp = SG.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
        
    %[V,Vdiff,VInt]   = vext(SG.Pts.y1_kv,SG.Pts.y2_kv);   
    [V,Vdiff,VInt]   = vext(SG.GetCartPts);   
    [VP]             = vext(SG.GetCartPts(Interp.pts1,Interp.pts2));
           
    %Check Differentiation
    %vplot   = Interp.InterPol*V;        
    %data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Interpolation
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);    
    
    %Check Integration
    data.Int      = abs(Int*V-VInt);
    data.IntValue = Int*V;
    display([' Error in Integration: ', num2str(data.Int)]);                                
    
%    subplot(2,2,1); doPlots_SC(Interp,Pts,Vdiff.dy1);
%    subplot(2,2,2); doPlots_SC(Interp,Pts,Diff.Dy1*V);
%    subplot(2,2,3); doPlots_SC(Interp,Pts,Diff.Dy1*V-Vdiff.dy1);    
               
    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    VI = V.*sin(Pts.y1_kv).*sin(Pts.y2_kv).^2*R^2;
    if(volume)
        VI = VI.*(2*R*sin(Pts.y1_kv).*sin(Pts.y2_kv));
    end
    SG.plot(VI,{'comp','SC'});
    SaveFigure(['SimulationBall_Sphere_Vol',num2str(volume)],shape);
    
    figure
    set(gcf,'Color','white'); %Set background color                
     SG.plot(V,'SC');        
    title('Interpolation');   
    
    figure;
    SG.PlotGridLines();    
    SG.PlotGrid();
    	
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;
    
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [V,VDiff,VInt] = VTest2(y1,y2)       
            %V           = cos(y1);  VInt = 0;
           % V           = sin(y1).^2.*cos(y2).^2;  VInt = 4/3*pi*R^2;
            V   = ones(size(y1));  VInt = 4/3*pi*R^3; 
           % V = 2*R*sin(y1).*sin(y2);
            
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

    function [V,VDiff,VInt] = VTest(pts)
        y1 = pts.y1_kv;
        y2 = pts.y2_kv;
        d  = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(1/2);        
        V  = BarkerHenderson_2D(d);        
        VDiff = 0;
        VInt = -32/9*pi+25/32*pi^2;
    end

    function [V,VDiff,VInt] = VTestt(pts)
        y1 = pts.y1_kv;
        y2 = pts.y2_kv;
        d  = ((y1-Origin(1)).^2+(y2-Origin(2)).^2).^(1/2);                
        rc = 5;
        
        V2.r_cutoff = rc;
        V2.epsilon  = 1;
        
        V  = BarkerHendersonCutoff_2D(d,V2);
        %V  = BarkerHenderson_2D(d);        
        VDiff = 0;
        if(rc == inf)
            VInt = -32/9*pi+25/32*pi^2;
        else
            VInt = pi * (0.18000e5 * atan(sqrt((rc ^ 2 - 1))) * (rc ^ 10) - (40960 * rc ^ 10) + (18000 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 6) - (9240 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 4) - 0.22200e5 * sqrt((rc ^ 2 - 1)) * (rc ^ 6) + (61440 * rc ^ 7) - (12978 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1) * rc ^ 2) + 0.1050e4 * sqrt((rc ^ 2 - 1)) * (rc ^ 4) - (16857 * (rc ^ 2 - 1) ^ (0.3e1 / 0.2e1)) + 0.1575e4 * sqrt((rc ^ 2 - 1)) * (rc ^ 2) + 0.1575e4 * sqrt((rc ^ 2 - 1)) - (20480 * rc)) / (rc ^ 10) / 0.11520e5;        
        end
    end

end