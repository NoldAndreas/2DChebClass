function data = SimulationBox(N1,N2,L1,L2,vext)

    disp('** Simulation Box **');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  15;    N2 = 15;
        L1    = 2;  L2    = 2;
        Origin = [-L1;-L2]/2;
        vext  = @Vext2;
        N     = [N1,N2];
    end        
    
    abox               = Box(v2struct(L1,L2,Origin,N));
    [Pts,Diff,Int,Ind] = abox.ComputeAll();
    
    Interp = abox.ComputeInterpolationMatrix((-1:0.02:1)',(-1:0.02:1)',true,true);
    Conv   = abox.ComputeConvolutionMatrix(@f2);
    
     %Check Spectral/Spectral map in 2D:    
    %[Pts,Diff,Int,Ind,Interp,Conv] = SpectralSpectral(Maps,N1,N2,y1Plot,y2Plot,@f2);                
    intBound = struct('y1_l',abox.PhysSpace1(-1),...
                      'y1_u',abox.PhysSpace1(1),...
                      'y2_l',abox.PhysSpace2(-1),...
                      'y2_u',abox.PhysSpace2(1));                               

    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound,'cart');    
    [VP]             = vext(Interp.pts1,Interp.pts2);                          

    %Check Differentiation
    vplot   = Interp.InterPol*V;
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
        
    %Check Convolution
    fP = f1(Pts.y1_kv,Pts.y2_kv);
    data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    display([' Error in Convolution: ', num2str(data.Conv)]);
    
    data.N1 = N1; data.N2 = N2;
%    data.Interp = Interp; data.f = V;    
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    abox.plot(V,'SC');     title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);
    abox.plot(fConv(Interp.pts1,Interp.pts2));
    title('Convolution');
    pbaspect([1 1 1]);

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(x,y)
        z = x;
    end
    function z = f2(r)%x,y)
        z = exp(-r.^2);%(x.^2+y.^2));
    end
    function z = fConv(x,y)
        z = (-sqrt(pi)*exp(-(x+1).^2) - pi*x.*erf(1+x) +...
              sqrt(pi)*exp(-(x-1).^2) + pi*x.*erf(-1+x)).*...
              (-erf(1+y)+erf(y-1))/4;
    end
end