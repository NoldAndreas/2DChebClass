function SimulationInfinityWedge

    disp('** SimulationPolarInfinityWedge **');
    
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    LR                   = 2;       
    vext                = @Vext9;
    N1                  = 30;
    N2                  = 20;
    
    R                   = 2;
    half_wedge_angle    = pi/4;
    Geometry = struct('R_in',0,'LR',LR,...
                      'th1',-half_wedge_angle,...
                      'th2',half_wedge_angle,...
                      'N',[N1,N2]);    
    
    tic
    WDG                = InfWedge(Geometry);
    x1plotMax          = WDG.CompSpace1(5);
    [Pts,Diff,Int,Ind] = WDG.ComputeAll();
    Interp = WDG.ComputeInterpolationMatrix((-1:0.02:x1plotMax)',...
                                                (-1:0.02:1)',true,true);
    Conv   = WDG.ComputeConvolutionMatrix(@f2,true);            
    toc       
    
    %Check Polar Spectral/Spectral map in 2D:   
         
    intBound = struct('y1_l',WDG.PhysSpace1(-1),...
                                  'y1_u',WDG.PhysSpace1(1),...
                                  'y2_l',WDG.PhysSpace2(-1),...
                                  'y2_u',WDG.PhysSpace2(1));
                              
    [V,Vdiff,VInt] = vext(Pts.y1_kv,Pts.y2_kv,intBound,'polar');    
    [VP]           = vext(Interp.pts1,Interp.pts2);               
    
    %Check Differentiation 
    vplot       = Interp.InterPol*V;        
    displayErrors(vplot,VP,V,Vdiff,Diff);
        
    %Check Interpolation   
    subplot(2,1,1);  WDG.doPlots(V);
    subplot(2,1,2);  WDG.doPlots(Interp.InterPol*V - VP);        

    %Check Integration
    display([' Error in Integration: ', num2str(Int*V-VInt)]);                    
    
    %Check Convolution
    figure
    fP_Conv  = Conv*f1(Pts.y1_kv,Pts.y2_kv);
    WDG.doPlots(fP_Conv);    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(r,t)
        x = r.*cos(t);
        y = r.*sin(t);
        
        x0 = 1;
        rn = sqrt((x-x0).^2+y.^2);
        
        s = 1;
        z = 1/(s*sqrt(2*pi))*exp(-(rn.^2)/(2*s^2));  
        z(r == inf) = 0;
    end
    function z = f2(r,t)       
        s = 2;
        z = 1/(s*sqrt(2*pi))*exp(-(r.^2)/(2*s^2));
        z(r == inf) = 0;        
    end
    

end

