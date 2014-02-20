function testInterpolationModel4
    close all;
    
    Maps    = struct('PhysSpace',@Comp_to_Phys,'CompSpace',@Phys_to_Comp);                
    Maps2   = struct('PhysSpace',@Comp_to_Phys2,'CompSpace',@Phys_to_Comp2);                
    yPlot   = (-100:0.02:100)';
    L       = 1;
    
    [h1,h2,fInt] =  f(0);
        
    n = (10:2:20)';
    l = (0.1:1:10)';
    
    N = 20;
    L = 1;
    for i = 1:length(n)        
    %    for i = 1:length(l)                        
            [errInterp(i), errDiff(i), errInt(i)]   = testInterpolation(Maps,n(i));            
            [errInterp2(i), errDiff2(i), errInt2(i)] = testInterpolation(Maps2,n(i));
    %    end
    end
    
    subplot(2,1,1);
    semilogy(n,errInterp,'b'); hold on; semilogy(n,errInterp2,'m'); 
    title('Error Interpolation');
    
    subplot(2,1,2);
    semilogy(n,errInt,'b'); hold on; semilogy(n,errInt2,'m'); 
    title('Error Integration');
    
    
    
    L = 10;
    [Pts,Diff,Int,Interp]     = Spectral(Maps,20,yPlot);        
    L =1;
    [Pts2,Diff2,Int2,Interp2] = Spectral(Maps2,20,yPlot);        
    [z,dz]     = f(Pts.y_kv);     zPlot   = f(yPlot);
    [z2,dz2]   = f(Pts2.y_kv);   z2Plot   = f(yPlot);
    
    %semilogy(Interp.pts,Interp.InterPol*z,'r'); semilogy(Pts.y_kv,z,'or'); hold on;
    %semilogy(Interp2.pts,Interp2.InterPol*z2,'g'); semilogy(Pts2.y_kv,z2,'og'); hold on;
    %semilogy(Interp.pts,zPlot,'b'); 
    
    semilogy(Interp.pts,abs(Interp.InterPol*z-zPlot),'r');% semilogy(Pts.y_kv,z,'or');
    hold on;
    semilogy(Interp2.pts,abs(Interp2.InterPol*z2-zPlot),'g'); %semilogy(Pts2.y_kv,z2,'og'); hold on;
    %semilogy(Interp.pts,zPlot,'b'); 
    
    
    
    
%     N = 15;
%     for j = 1:length(l)
%         L = l(j);
%         [errInterp(j), errInt(j)]   = testInterpolation(Maps,N);
%         [errInterp2(j), errInt2(j)] = testInterpolation(Maps2,N);
%     end
%     
%     plot(l,errInterp); hold on;
%     plot(l,errInterp2,'m'); hold on;
%     
    
    
    for i = 1:length(n)
        disp(num2str(n(i)));
        for j = 1:length(l)
            L = l(j);
            [errInterp(i,j), errDiff(i,j), errInt(i,j)]   = testInterpolation(Maps,n(i));
            [errInterp2(i,j), errDiff2(i,j), errInt2(i,j)] = testInterpolation(Maps2,n(i));
        end
    end
    nMat = repmat(n,1,length(l));
    lMat = repmat(l',length(n),1); 
    
    
    subplot(3,2,1);
    contourf(nMat,lMat,log(errInterp)/log(10));
    cMin = min(min([log(errInterp)/log(10);log(errInterp2)/log(10)]));
    cMax = max(max([log(errInterp)/log(10);log(errInterp2)/log(10)]));
    caxis([cMin cMax]);   
    xlabel('N'); ylabel('L');
    title('Interpolation');    
    
    subplot(3,2,2);
    contourf(nMat,lMat,log(errInterp2)/log(10));
    xlabel('N'); ylabel('L'); title('Interpolation');
    colorbar; 
    cMin = min(min([log(errInterp)/log(10);log(errInterp2)/log(10)]));
    cMax = max(max([log(errInterp)/log(10);log(errInterp2)/log(10)]));
    caxis([cMin cMax]);        
    
    subplot(3,2,3);
    contourf(nMat,lMat,log(errInt)/log(10));
    cMin = min(min([log(errInt)/log(10);log(errInt2)/log(10)]));
    cMax = max(max([log(errInt)/log(10);log(errInt2)/log(10)]));
    colorbar;  caxis([cMin cMax]);           
    xlabel('N'); ylabel('L'); title('Integration');     
    
    subplot(3,2,4);
    contourf(nMat,lMat,log(errInt2)/log(10));
    xlabel('N'); ylabel('L');  title('Integration');
    cMin = min(min([log(errInt)/log(10);log(errInt2)/log(10)]));
    cMax = max(max([log(errInt)/log(10);log(errInt2)/log(10)]));
    colorbar;  caxis([cMin cMax]);        
    
    subplot(3,2,5);
    contourf(nMat,lMat,log(errDiff)/log(10));
    xlabel('N'); ylabel('L'); title('Differentiation');
    cMin = min(min([log(errDiff)/log(10);log(errDiff2)/log(10)]));
    cMax = max(max([log(errDiff)/log(10);log(errDiff2)/log(10)]));
    colorbar;  caxis([cMin cMax]);        
        
    subplot(3,2,6);
    contourf(nMat,lMat,log(errDiff2)/log(10));
    xlabel('N'); ylabel('L'); title('Differentiation');
    cMin = min(min([log(errDiff)/log(10);log(errDiff2)/log(10)]));
    cMax = max(max([log(errDiff)/log(10);log(errDiff2)/log(10)]));
    colorbar;  caxis([cMin cMax]);        
    
    subplot(2,1,1);
    semilogy(repmat(n,1,length(l)),errInterp);
    for j = 1:length(l)        
        if j == 1
            str = {num2str(l(j))};
        else
            str = [str, num2str(l(j))];
        end
    end
    legend(str);
    
    %hold on; semilogy(n,errInterp2,'m'); 
    title('Error Interpolation');        
    
    subplot(2,1,2);
    semilogy(n,errIn t,'b'); hold on; semilogy(n,errInt2,'m'); 
    title('Error Integration');
    
    function [errInterp,errDiff, errInt] = testInterpolation(MapsV,N)
        
        [Pts,Diff,Int,Interp] = Spectral(MapsV,N,yPlot);
        
        [zt,dzt]   = f(Pts.y_kv);   zPlot   = f(yPlot);
        
        %semilogy(Interp.pts,Interp.InterPol*z); hold on;
        %semilogy(Interp.pts,zPlot,'k'); hold off;
        %doPlots_SC_1D(Interp,Pts,z); hold on; doPlots_IP_1D(Interp,zPlot);
        errInterp = max(abs(Interp.InterPol*zt - zPlot));
        errInt    = abs(Int*zt - fInt);
        errDiff   = max(abs(dzt - Diff.Dy*zt));
        %disp(['Interp/Integr Error: ', num2str(errInterp),' , ', num2str(errInt)]);         
        
    end
    
    %****************************************
    function [z,dz,dx,ddx,dddx] = Comp_to_Phys(x)
        [z,dz,dx,ddx,dddx] = SqrtMap(x,L,inf); 
    end
    function xf = Phys_to_Comp(z)
        xf = InvSqrtMap(z,L,inf); 
    end


    function [z,dz,dx,ddx,dddx] = Comp_to_Phys2(x)
        [z,dz,dx,ddx,dddx] = MapInf2(x,L/sqrt(2)); 
    end
    function xf = Phys_to_Comp2(z)
        xf = InvMapInf2(z,L/sqrt(2)); 
    end

    function [y,dy,yInt] = f(r)
        a = 1;
        %y       = 1./((a^2+r.^2)); yInt    = pi/a;  dy = -2*y./((a^2+r.^2).^2);
       % y       = 1./((a^2+r.^2).^(3/2));  dy      = -3*r./((a^2+r.^2).^(5/2));  yInt    = 2/a^2;
        y       = 1./((1+r.^2).^2);     yInt    = pi/2; dy = -4*y./((1+r.^2).^3);
        % y       = 1./((1+r.^2).^(5/2)); yInt    = 4/3;
        %y = exp(-r.^2); yInt = sqrt(pi);
    end

end