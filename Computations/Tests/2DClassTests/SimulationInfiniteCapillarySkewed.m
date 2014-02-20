function data = SimulationInfiniteCapillarySkewed()

    disp('** Simulation Infinite Capillary Skewed **');

    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %**************************************************
    %Parameters
    N1 = 30;  L1    = 2; 
    N2 = 40;  L2    = 1;
    N  = [N1,N2];
    alpha = pi/3;

    Conv   = struct('L1',L1,'L2',L2,'N',[N1,N2]);
    VextC  = @Vext7;
    y2Min  = 0; 
    y2Max  = 1;
    %**************************************************
    
    
    IC     = InfCapillarySkewed(v2struct(L1,y2Min,y2Max,N,Conv,alpha));
    x1PAbs = IC.CompSpace1(3);
    [Pts,Diff,Int,Ind] = IC.ComputeAll();    
    Interp = IC.ComputeInterpolationMatrix((-x1PAbs:0.02:x1PAbs)',(-1:0.02:1)',true,true);
    Conv   = IC.ComputeConvolutionMatrix(@f2,true);                         
    
    intBound = struct('y1_l',IC.PhysSpace1(-1),...
                      'y1_u',IC.PhysSpace1(1),...
                      'y2_l',IC.PhysSpace2(-1),...
                      'y2_u',IC.PhysSpace2(1));
                       
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound);    
    [VP]             = vext(Interp.pts1,Interp.pts2);                          
    
    %Check Differentiation
    vplot     = Interp.InterPol*V;        
    data      = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                

    data.N1 = N1; data.N2 = N2;
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(2,1,1);
    IC.doPlots(V,'SC');    
    title('Interpolation');    
    %pbaspect([1 1 1]);
    
    subplot(2,1,2);
    IC.doPlots(Conv(:,Pts.y1_kv < inf)*V(Pts.y1_kv < inf),'SC');
    title('Convolution');
    %pbaspect([1 1 1]);
                                    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f2(r)%x,y)
        z = exp(-r.^2);%(x.^2+y.^2));
    end

    function [V,Vdiff,VInt] = vext(y1t,y2t,intBound)
        ptsCart        = IC.GetCartPts(y1t,y2t);        
        
        if(nargin>2)
            [V,Vdiff,VInt] = VextC(ptsCart.y1_kv,ptsCart.y2_kv,intBound,'cart');
        else
            [V,Vdiff,VInt] = VextC(ptsCart.y1_kv,ptsCart.y2_kv);            
        end
    end
end