function data = SimulationPolarInfinity(N1,N2,L1,L2,vext)

    disp('** SimulationPolarInfinityClass **');
    
%    if(length(dbstack) == 1)
%        AddPaths();
%    end    
    close all;
    
    %Initialization    
    if(nargin == 0)
        L       = 2;
        vext    = @Vext7; %@Vext10;
        N1      = 20;
        N2      = 20;
        N       = [N1;N2];
    else
        L = L1;
    end

    IDC                = InfDisc(v2struct(L,N));
    x1plotMax          = IDC.CompSpace1(10);                      
    [Pts,Diff,Int,Ind] = IDC.ComputeAll();    
    Interp = IDC.ComputeInterpolationMatrix((0:100)'*x1plotMax/100,...
                                                (0:0.02:1)',true,true);
                                            
    disp('The convolution test has been deactivated, because now only functions which depend on the distance are accepted.');                                            
    %IDC.ComputeConvolutionMatrix_Test();    
                                            
    shapeParams = struct('L',3,'N',[10;10]);
    Conv   = IDC.ComputeConvolutionMatrix(@f2,shapeParams);
    
    intBound         = struct('y1_l',IDC.PhysSpace1(0),...
                              'y1_u',IDC.PhysSpace1(1),...
                              'y2_l',IDC.PhysSpace2(0),...
                              'y2_u',IDC.PhysSpace2(1));
    
    [V,Vdiff,VInt] = vext(Pts.y1_kv,Pts.y2_kv,intBound,'polar');    
    [VP]           = vext(Interp.pts1,Interp.pts2); 
    
    %Test Differentiation
    vplot   = real(Interp.InterPol*V);
    data    = displayErrors(vplot,VP,V,Vdiff,Diff);	                         
                
    %Test integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);
    
    %Test Interpolation
    data.InterPol = max(abs(Interp.InterPol*V - VP));       
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
                    
    %Test Convolution    
    fP_Conv   = Conv*f1(Pts.y1_kv,Pts.y2_kv);
    data.Conv = max(abs(Interp.InterPol*fP_Conv-fConv(Interp.pts1,Interp.pts2)));                     
    display([' Error in Convolution: ', num2str(data.Conv)]);    
    
    data.N1 = N1; data.N2 = N2;  
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    IDC.doPlots(V); 
    title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);    
    IDC.doPlots(fP_Conv);
    title('Convolution');
    pbaspect([1 1 1]);            
    
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
    function z = fConv(r,t)
        x  = r.*cos(t);
        y  = r.*sin(t);
        
        s  = sqrt(1+4);
        x0 = 1;
        rn = sqrt((x-x0).^2+y.^2);
                
        z = 2/(s^2)*exp(-(rn.^2)/(2*s^2));
        z(r== inf) = 0;
    end

end

