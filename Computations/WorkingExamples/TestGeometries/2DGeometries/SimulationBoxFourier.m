function [PB,res] = SimulationBoxFourier

    disp('** Simulation Box Fourier **');
    
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    L1      = 2;
    L2      = 2*pi;
    Origin  = [-L1/2,0];    %[0;0];
    N       = [20;20];    
    
    PB                 = PeriodicBox(v2struct(L1,L2,Origin,N));
    [Pts,Diff,Int,Ind] = PB.ComputeAll();
    
    Interp = PB.ComputeInterpolationMatrix((-1:0.02:1)',(0:0.02:1)',true,true);
    Conv   = PB.ComputeConvolutionMatrix(@f2);        
    
     %Check Spectral/Fourier map in 2D:
     %[Pts,Diff,Int,Ind,Interp,Conv] = SpectralFourier(Maps,n1,n2,rPlot,tPlot,@f2);
         
     intBound = struct('y1_l',PB.PhysSpace1(-1),...
                       'y1_u',PB.PhysSpace1(1),...
                       'y2_l',PB.PhysSpace2(0),...
                       'y2_u',PB.PhysSpace2(1));                               

     [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound);    
     [VP]             = vext(Interp.pts1,Interp.pts2);           
            
     %Check Differentiation and Interpolation
     displayErrorsPos(Pts,Interp.InterPol*V,VP,V,Vdiff,Diff,'cart');        
                 
     subplot(2,1,1); PB.plot(V);
     subplot(2,1,2); PB.plot(Interp.InterPol*V-VP);    
     title('Error');
     set(gcf,'name','Interpolation');
     
     %Check Integration
     display([' Error in Integration: ', num2str(Int*V-VInt)]);                
    
     %Check Convolution               
     figure
     fC_fft     = Conv*f1(Pts.y1_kv,Pts.y2_kv);
     subplot(2,1,1); PB.plot(fC_fft);
     subplot(2,1,2); PB.plot(Interp.InterPol*fC_fft-fConv(Interp.pts1,Interp.pts2));          
     title('Error');
     set(gcf,'name','Convolution');                                                                
     
     figure;
    PB.PlotGridLines();    
    PB.PlotGrid();
    	
    hl = xlabel('$y_1$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);
    hl = ylabel('$y_2$'); set(hl,'Interpreter','Latex'); set(hl,'fontsize',25);        
    res.fig_handles{1} = gcf;

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(x,y)
        %z = (x.^2).*cos(y);
        z = x.*sin(y);
    end
    function z = f2(x,y)
        z = cos(y).*exp(-x.^2);
        %z = x.*(1+sin(y));
    end
    function z = fConv(x,y)
        %z = 2/3*pi*x.*sin(y);  
        if(Origin(1) == -1)
            z = pi*sin(y)/2.*...
                (exp(-(1+x).^2) - exp(-(1-x).^2) +...
                sqrt(pi)*(x.*erf(1+x) - x.*erf(x-1)));
        elseif(Origin(1) == 0)
            z = pi*sin(y)/2.*...
                (x.*sqrt(pi).*(erf(x)-erf(x-2)) + exp(-x.^2) - exp(-(x-2).^2));
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

            k = 4;
            V           = 1+ (y1.^2).*(sin(k*y2));  
            dVr         = 2*y1.*sin(k*y2);        
            ddVr        = 2*sin(k*y2);        
            dVt         = k*(y1.^2).*cos(k*y2);  
            ddVt        = -k^2*(y1.^2).*sin(k*y2);  
            dVrdVt      = 2*k*y1.*cos(k*y2);       

            VDiff       = struct('dy1',dVr,'dy2',dVt,'ddy1',ddVr,'ddy2',ddVt,...
                                'dy1dy2',dVrdVt);

            VInt        = 0;
            if(nargin > 2)
                vr = vextIntR(intBound.y1_u) - vextIntR(intBound.y1_l);
                vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
                VInt = vr*vt + (intBound.y1_u - intBound.y1_l)*(intBound.y2_u - intBound.y2_l);
            end
        function vi = vextIntR(r)
                vi = r.^3/3;       
        end

        function vi = vextIntT(t)
            vi = cos(k*t)/k;
        end

    end
end