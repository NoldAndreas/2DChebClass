function testBoyd

    close all;
    MapsRatMap    = struct('PhysSpace',@Comp_to_Phys,'CompSpace',@Phys_to_Comp);                    
    MapsQuadMap   = struct('PhysSpace',@Comp_to_Phys2,'CompSpace',@Phys_to_Comp2);                
    yPlot         = (-10:0.1:10)';
    
    L = 4;
    
    
    n = 10:2:40;
    row  = 1;
    rows = 4;
    
    figure;
    set(gcf,'Color','white'); %Set background color    
    ff1 = figure(1);
    set(ff1, 'Position', [0 0 850 1500]);%3000 2000] );    
    %set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    
    
    
    PlotErrors(@f2,1,'e^{-y^2}');
    PlotErrors(@f3,2,'(1+y^2)^{-2}');       
    PlotErrors(@f1,3,'(1+y^2)^{-3/2}');
    PlotErrors(@f4,4,'(1+y^2)^{-5/2}');
    
    function PlotErrors(ft,row,str)
        f = ft;
    for j = 1:length(n)
        [errInterpQuadMap(j), errDiffQuadMap(j), errIntQuadMap(j)] = getErrorQuadMap(n(j),ft);
        [errInterpRatMap(j), errDiffRatMap(j), errIntRatMap(j)]    = getErrorRatMap(n(j),ft);
    end
    
    
    subplot(rows,3,(row-1)*3 + 1);
    semilogy(n,errIntQuadMap,'b'); hold on;
    semilogy(n,errIntRatMap,'m--');
    legend('QuadMap','RatMap');
    title('Error of Integration');
    pbaspect([1 1 1]);
    xlabel('N');
    %ylabel(str,'Interpreter','Latex');
    ylabel(str);
    
    subplot(rows,3,(row-1)*3 +2);
    semilogy(n,errInterpQuadMap,'b'); hold on;
    semilogy(n,errInterpRatMap,'m--');
    legend('QuadMap','RatMap');
    title('Max. Error of Interpolation');
    pbaspect([1 1 1]);
    
    subplot(rows,3,(row-1)*3 +3);
    semilogy(n,errDiffQuadMap,'b'); hold on;
    semilogy(n,errDiffRatMap,'m--');
    legend('QuadMap','RatMap');
    title('Max. Error of Differentiation');
    pbaspect([1 1 1]);
    
    end
    
    function [errInterpol,errDiff,errInt] = getErrorQuadMap(N,f)
    
        %i        = 0:(N+1);
        %ti       = i*pi/(N+1);
        %yi       = L*cot(ti);
        %wi       = [0 (L*pi./((N+1)*(sin(ti(2:end-1))).^2)) 0];
        %[z,zInt] = f(yi');
        %IntV     = wi*z;% + pi/(L*(N+1));
        
        [Pts,Diff,Int,Interp] = Spectral(MapsQuadMap,N,yPlot);
        [z,dz,ddz,zInt] = f(Pts.y_kv);
        zInterpol       = f(yPlot);
        IntV = Int*z;                
        
        errInt      = abs(IntV- zInt);
        errInterpol = max(abs(zInterpol-Interp.InterPol*z));
        %errDiff     = max(abs(dz-Diff.Dy*z));
        errDiff     = max(abs(ddz-Diff.DDy*z));        
    end

    function [errInterpol,errDiff,errInt] = getErrorRatMap(N,f)
    
        %i        = 0:(N+1);
        %ti       = i*pi/(N+1);
        %yi       = L*cot(ti);
        %wi       = [0 (L*pi./((N+1)*(sin(ti(2:end-1))).^2)) 0];
        %[z,zInt] = f(yi');
        %IntV     = wi*z;% + pi/(L*(N+1));
        
        [Pts,Diff,Int,Interp] = Spectral(MapsRatMap,N,yPlot);
        i               = 0:(N-1);
        ti              = i*pi/(N-1);        
        Int             = [0 (L*pi./((N-1)*(sin(ti(2:end-1))).^2)) 0];        
        %Int            = [0 (L*pi./((N+1)*(sin(ti(2:end-1))).^2)) 0];
        [z,dz,ddz,zInt] = f(Pts.y_kv);
        zInterpol       = f(yPlot);
        IntV            = Int*z;                
        
        errInt      = max(abs(IntV- zInt),eps);
        errInterpol = max(max(abs(zInterpol-Interp.InterPol*z)),eps);
        %errDiff     = max(abs(dz-Diff.Dy*z));        
        errDiff     = max(abs(ddz-Diff.DDy*z));        

    end

    %function [z,zInt] = f(y)
        %z = 1./(1+(y+1).^2); zInt = pi;
 %       z = 1./((1+y.^2).^(3/2)); zInt = 2;
        %z = 1./((1+y.^2).^2); zInt = pi/2;
        %z = 1./((1+y.^2).^(5/2)); zInt = 4/3;
        %z = 1./(1+y.^2).^2; zInt = pi/2;
        %z  = exp(-y.^2); zInt = sqrt(pi);
%    end

    function [y,dy,ddy,yInt] = f1(r)
        a    = 1;
        y    = 1./((a^2+r.^2).^(3/2));  
        dy   = -3*r./((a^2+r.^2).^(5/2));  
        ddy  = 3*(4*r.^2-a^2)./((a^2+r.^2).^(7/2));  
        yInt = 2/a^2;
    end
    function [y,dy,ddy,yInt] = f2(r)                       
        y    = exp(-r.^2); 
        dy   = -2*r.*y; 
        ddy  = (-2+4*r.^2).*y; 
        yInt = sqrt(pi);
    end
    function [y,dy,ddy,yInt] = f3(r)
       y    = 1./((1+r.^2).^2); 
       yInt = pi/2; 
       dy   = -4*r./((1+r.^2).^3);         
       ddy  = 4*(5*r.^2-1)./((1+r.^2).^4);         
    end
    function [y,dy,ddy,yInt] = f4(r)
       y    = 1./((1+r.^2).^(5/2)); 
       yInt = 4/3; 
       dy   = -5*r./((1+r.^2).^(7/2));        
       ddy  = 5*(6*r.^2-1)./((1+r.^2).^(9/2));        
    end


    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys(x)        
        [z,dz,dx,ddx,dddx,ddddx] = SqrtMap(x,L,inf); 
        %z   = L*x./((1-x.^2).^(1/2));
        %dz  = L./((1-x.^2).^(3/2));
    end
    function xf = Phys_to_Comp(z)
        xf = InvSqrtMap(z,L,inf); 
    end

    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys2(x)
        [z,dz,dx,ddx,dddx,ddddx] = QuadMap(x,L,inf); 
        %[z,dz,dx,ddx,dddx,ddddx] = QuadMapAB(x,L,-5,5); 
    end
    function xf = Phys_to_Comp2(z)
        xf = InvQuadMap(z,L,inf); 
        %xf = InvQuadMapAB(z,L,-5,5); 
    end


end