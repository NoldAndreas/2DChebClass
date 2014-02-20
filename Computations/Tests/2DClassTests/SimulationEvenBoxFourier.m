function SimulationEvenBoxFourier

    disp('** Simulation Even Box Fourier **');
    
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    L1      = 2;
    L2      = 2*pi;
    Origin  = [-1;0];    
    N       = [20;20];    
    
    PB                 = EvenPeriodicBox(v2struct(L1,L2,Origin,N));
    [Pts,Diff,Int,Ind] = PB.ComputeAll();
    
    Interp = PB.ComputeInterpolationMatrix((-1:0.02:1)',(0:0.02:1)',true,true);
    Conv   = PB.ComputeConvolutionMatrix(@f2);        
    
     %Check Spectral/Fourier map in 2D:
     %[Pts,Diff,Int,Ind,Interp,Conv] = SpectralFourier(Maps,n1,n2,rPlot,tPlot,@f2);
         
     intBound = struct('y1_l',PB.PhysSpace1(-1),...
                       'y1_u',PB.PhysSpace1(1),...
                       'y2_l',PB.PhysSpace2(0),...
                       'y2_u',PB.PhysSpace2(0.5));                               

     [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound);    
     [VP]             = vext(Interp.pts1,Interp.pts2);           
            
     %Check Differentiation and Interpolation
     displayErrorsPos(Pts,Interp.InterPol*V,VP,V,Vdiff,Diff,'cart');        
                 
     subplot(2,1,1); PB.doPlots(V);
     subplot(2,1,2); PB.doPlots(Interp.InterPol*V-VP);    
     title('Error');
     set(gcf,'name','Interpolation');
     
     %Check Integration
     display([' Error in Integration: ', num2str(Int*V-VInt)]);                
    
     %Check Convolution               
     figure
     fC_fft     = Conv*f1(Pts.y1_kv,Pts.y2_kv);
     subplot(2,1,1); PB.doPlots(fC_fft);
     subplot(2,1,2); PB.doPlots(Interp.InterPol*fC_fft-fConv(Interp.pts1,Interp.pts2));          
     title('Error');
     set(gcf,'name','Convolution');                                                                

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(x,y)
        %z = (x.^2).*cos(y);
        z = x.*cos(y);
    end
    function z = f2(x,y)
        z = cos(y).*exp(-x.^2);
        %z = x.*(1+sin(y));
    end
    function z = fConv(x,y)
        %z = 2/3*pi*x.*sin(y);        
        if(Origin(1) == -1)
            z = pi*cos(y)/2.*...
                (exp(-(1+x).^2) - exp(-(1-x).^2) +...
                sqrt(pi)*(x.*erf(1+x) - x.*erf(x-1)))/2;
        elseif(Origin(1) == 0)
            z = pi*cos(y)/2.*...
                (x.*sqrt(pi).*(erf(x)-erf(x-2)) + exp(-x.^2) - exp(-(x-2).^2))/2;
        else
            disp('**************Error*******************');
            disp('Convolution will not be computed correctly');
        end        
    end


    function [V,VDiff,VInt] = vext(y1,y2,intBound)       

        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u       
        
        V           = 1+ (y1.^2).*cos(y2);  
        dVr         = 2*y1.*cos(y2);        
        dVt         = -(y1.^2).*sin(y2);  
        
        VDiff       = struct('dy1',dVr,'dy2',dVt);
        
        VInt = 0;
        if(nargin > 2)
            vr = vextIntR(intBound.y1_u) - vextIntR(intBound.y1_l);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vr*vt + (intBound.y1_u - intBound.y1_l)*(intBound.y2_u - intBound.y2_l);
        end
        
        function vi = vextIntR(r)
            vi = (r.^4)/4;
        end

        function vi = vextIntT(t)
            vi = sin(t);
        end

    end

end