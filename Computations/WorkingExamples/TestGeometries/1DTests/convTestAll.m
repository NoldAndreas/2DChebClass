function convTestAll

    t1 = 5;
    t2 = 25;

    shape.L = 5;
    shape.N = 20;
    
%    plotRange = struct('yMin',-5,'yMax',5,'N',100);
    
    aline               = InfSpectralLine(shape);
%    [Pts,Diff,Int,Ind,Interp] = aline.ComputeAll(plotRange);    

    y = aline.Pts.y;
    g = Gaussian(y);
    
    function g=Gaussian(x)
        g = exp(-x.^2);
    end


    function k = kernel(x,t)
        k = exp(-x.^2/4/t);
    end

    function k = kernels(x)
        k = [kernel(x,t1) kernel(x,t2)];
    end

    
    %----------------------------------------------------------------------
    % Full space
    %----------------------------------------------------------------------
    
    infGeom.N = 100; infGeom.L = 5;
    M_Inf = aline.ComputeConvolutionMatrix(@kernels,infGeom,false);

    figure
    cInf1 = M_Inf(:,:,1)*g;
    cInf2 = M_Inf(:,:,2)*g;
    plot(y,cInf1,'b')
    hold on
    plot(y,cInf2,'r')
    
    plot(y,exactInf(y,t1),'ob');
    plot(y,exactInf(y,t2),'or'); 
    
    function f = exactInf(x,t)
        f = 2 * exp(-x.^2/(1+4*t)) * sqrt(pi/(4+1/t));
    end
    
    %----------------------------------------------------------------------
    % Finite space
    %----------------------------------------------------------------------
    
    yMax = 6;
    fGeom.N = 100; fGeom.yMax = yMax;
    M_f = aline.ComputeConvolutionMatrix(@kernels,fGeom,false);

    figure
    cf1 = M_f(:,:,1)*g;
    cf2 = M_f(:,:,2)*g;
    plot(y,cf1,'b')
    hold on
    plot(y,cf2,'r')
    
    plot(y,exactf(y,t1),'ob');
    plot(y,exactf(y,t2),'or'); 
    
    function f = exactf(x,t)
        f = exp(-x.^2/(1+4*t)) .* sqrt(pi*t/(1+4*t)) ...
                .* ( erf( (yMax - 4*t*(-yMax+x))/(2*sqrt(t*(1+4*t))) ) ...
                 - erf( (-yMax - 4*t*(yMax+x))/(2*sqrt(t*(1+4*t))) ) );
    end

    %----------------------------------------------------------------------
    % Infinite with gap
    %----------------------------------------------------------------------
    
    yMin = 3;
    
    infGapGeom.N = 100; infGapGeom.yMin= yMin; infGapGeom.L = 5;
   
    M_InfGap = aline.ComputeConvolutionMatrix(@kernels,infGapGeom,false);

    figure
    cInfGap1 = M_InfGap(:,:,1)*g;
    cInfGap2 = M_InfGap(:,:,2)*g;
    plot(y,cInfGap1,'b')
    hold on
    plot(y,cInfGap2,'r')
    
    plot(y,exactInfGap(y,t1),'ob');
    plot(y,exactInfGap(y,t2),'or'); 
    
    function f = exactInfGap(x,t)
        f_full =  2 * exp(-x.^2/(1+4*t)) * sqrt(pi/(4+1/t));
        f_exclude = exp(-x.^2/(1+4*t)) .* sqrt(pi*t/(1+4*t)) ...
             .* ( erf( (yMin - 4*t*(-yMin+x))/(2*sqrt(t*(1+4*t))) ) ...
                 - erf( (-yMin - 4*t*(yMin+x))/(2*sqrt(t*(1+4*t))) ) );
        f = f_full - f_exclude;
    end

    %----------------------------------------------------------------------
    % Finite with gap
    %----------------------------------------------------------------------
    
    yMin = 3;
    yMax = 6;
    
    fGapGeom.N = 100; fGapGeom.yMin= yMin; fGapGeom.yMax= yMax;
   
    M_fGap = aline.ComputeConvolutionMatrix(@kernels,fGapGeom,false);

    figure
    cfGap1 = M_fGap(:,:,1)*g;
    cfGap2 = M_fGap(:,:,2)*g;
    plot(y,cfGap1,'b')
    hold on
    plot(y,cfGap2,'r')
    
    plot(y,exactfGap(y,t1),'ob');
    plot(y,exactfGap(y,t2),'or'); 

    function f = exactfGap(x,t)
        f_out = exp(-x.^2/(1+4*t)) .* sqrt(pi*t/(1+4*t)) ...
             .* ( erf( (yMax - 4*t*(-yMax+x))/(2*sqrt(t*(1+4*t))) ) ...
                 - erf( (-yMax - 4*t*(yMax+x))/(2*sqrt(t*(1+4*t))) ) );
        f_in = exp(-x.^2/(1+4*t)) .* sqrt(pi*t/(1+4*t)) ...
             .* ( erf( (yMin - 4*t*(-yMin+x))/(2*sqrt(t*(1+4*t))) ) ...
                 - erf( (-yMin - 4*t*(yMin+x))/(2*sqrt(t*(1+4*t))) ) );
         f = f_out-f_in;
    end

   

end