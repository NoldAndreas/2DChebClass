function data = SimulationWedge_M1(N1,N2,R,vext)

    disp('** Simulation Wedge M1 **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  30;   N2 = 20; N = [N1;N2];
        R = 20;
        vext  = @Vext13;
        th1 = 0; th2 = pi;
    end                
    
    intBound = struct('y1_l',-1,'y1_u',1,'y2_l',-1,'y2_u',1);    
    
    %Check Spectral/Spectral map in 2D:
    WI     = Wedge_Interface(v2struct(R,th1,th2,N));
    [Pts,Diff,Int,Ind] = WI.ComputeAll();
    Interp             = WI.ComputeInterpolationMatrix((-1:0.05:1)',...
                                (-1:0.05:1)',true,true);                    
       
    
% %     Interp_r2   = M1SpectralSpectral_Interpolation(Pts.x1(2),x2Plot,Pts,Maps);     
% %     Interp_r4   = M1SpectralSpectral_Interpolation(Pts.x1(4),x2Plot,Pts,Maps);
% %     
% %     Interp_t0   = M1SpectralSpectral_Interpolation(x1Plot,0,Pts,Maps);     
% %     Interp_tPiD2   = M1SpectralSpectral_Interpolation(x1Plot,pi/2,Pts,Maps);     
% %     Interp_tPi   = M1SpectralSpectral_Interpolation(x1Plot,pi,Pts,Maps);     
    
%     TEST
%     PtsP.x1 = (-1:0.001:-0.9)';
%     PtsP.x2 = (-1:0.02:1)';
%     x1P_kv       = kron(PtsP.x1,ones(size(PtsP.x2)));
%     x2P_kv       = kron(ones(size(PtsP.x1)),PtsP.x2);
%     PtsP.N1 = length(PtsP.x1); PtsP.N2 = length(PtsP.x2);
%     [y1_kv,y2_kv,J,dH1,dH2] = Comp_to_Phys(x1P_kv,x2P_kv);
%     
%     subplot(2,2,1);     doPlots_Comp(PtsP,y2_kv);       title('y_2');
%     subplot(2,2,2);     doPlots_Comp(PtsP,J(:,2,1));    title('dy_2/dx_1');
%     subplot(2,2,3);     doPlots_Comp(PtsP,J(:,2,2));    title('dy_2/dx_2');
%     subplot(2,2,4);     doPlots_Comp(PtsP,dH2(:,1,1));  title('d^2y_2/dx_1^2');
%     
%     subplot(2,2,1);  doPlots_SC_Polar(Interp,Pts,dH2(:,1,1));
%     subplot(2,2,2);  doPlots_SC_Polar(Interp,Pts,dH2(:,1,2));
%     subplot(2,2,3);  doPlots_SC_Polar(Interp,Pts,dH2(:,2,1));
%     subplot(2,2,4);  doPlots_SC_Polar(Interp,Pts,dH2(:,2,2));
    
%    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound);    
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]             = vext(Interp.pts1,Interp.pts2);                          
    
    %doPlots_SC_Polar(Interp,Pts,V);
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
   
%    legend('r=0',['r=',num2str(Interp_r1.pts1(1))],['r=',num2str(Interp_r4.pts1(1))]);
    WI.plot(V);
    
    subplot(2,2,1); WI.plot(Vdiff.dy1);
    subplot(2,2,2); WI.plot(Diff.Dy1*V);
    subplot(2,2,3); WI.plot((Diff.Dy1*V-Vdiff.dy1));
    
    figure    
    
    subplot(2,2,1); WI.plot(Vdiff.ddy1);
    subplot(2,2,2); WI.plot(Diff.DDy1*V);
    subplot(2,2,3); WI.plot((Diff.DDy1*V-Vdiff.ddy1));
    
    figure

    subplot(2,2,1); WI.plot(Vdiff.dy1dy2);
    subplot(2,2,2); WI.plot(Diff.Dy1Dy2*V);
    subplot(2,2,3); WI.plot((Diff.Dy1Dy2*V-Vdiff.dy1dy2));

    rM0 = 1./Pts.y1_kv;
    rM0(rM0 == inf) = 0;   
    
    figure
        
    subplot(2,2,1); WI.plot(Vdiff.ddy2);
    subplot(2,2,2); WI.plot(Diff.DDy2*V);
    subplot(2,2,3); WI.plot((Diff.DDy2*V-Vdiff.ddy2));

    figure
    subplot(2,2,1); WI.plot(Vdiff.Lap);
    subplot(2,2,2); WI.plot(Diff.Lap*V);
    subplot(2,2,3); WI.plot(Diff.Lap*V-Vdiff.Lap);

    %Check Integration
    data.Int = abs(Int*V-VInt);
    display([' Error in Integration: ', num2str(data.Int)]);                
        
    %Check Convolution
    %fP = f1(Pts.y1_kv,Pts.y2_kv);
    %data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    %display([' Error in Convolution: ', num2str(data.Conv)]);
    
    data.N1 = N1; data.N2 = N2;
%    data.Interp = Interp; data.f = V;    
    
%    Plot1DComp_r(V,1);
%    Plot1DComp_r(V,4);
    Plot1DComp_r(V,11,0);
    Plot1DComp_r(V,11,1);
    Plot1DComp_r(V,11,2);
    
    figure
    Plot1DComp_t(V,1);
    Plot1DComp_t(V,4);
    Plot1DComp_t(V,10);

    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    
    WI.plot(V); 
    title('Interpolation');    
    pbaspect([1 1 1]);        

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f1(x,y)
        z = x;
    end
    function z = f2(x,y)
        z = exp(-(x.^2+y.^2));
    end
    function z = fConv(x,y)
        z = (-sqrt(pi)*exp(-(x+1).^2) - pi*x.*erf(1+x) +...
              sqrt(pi)*exp(-(x-1).^2) + pi*x.*erf(-1+x)).*...
              (-erf(1+y)+erf(y-1))/4;
    end

    function Plot1DComp_t(V,index)
        
        
        Interp_r   = M1SpectralSpectral_Interpolation(Pts.x1(index),x2Plot,Pts,Maps);     
        mark       = (Pts.y1_kv == Pts.y1_kv(N2*(index-1)+1));
        [~,Vdiffh,~]   = vext(Pts.y1_kv(mark),Pts.y2_kv(mark));        
        
        diffM  = Diff.DDy1;   
        VdiffD = Vdiffh.ddy1;
        
        plot(Interp_r.pts2,Interp_r.InterPol*(diffM*V),'-'); hold on;
        plot(Pts.y2_kv(mark),diffM(mark,:)*V,'o'); hold on;
        plot(Pts.y2_kv(mark),VdiffD,'-.');
    end
    function Plot1DComp_r(V,index,order)

        
        Interp_r   = M1SpectralSpectral_Interpolation(x1Plot,Pts.x2(index),Pts,Maps);     
        x2_kv     = kron(ones(size(Pts.x1)),Pts.x2);
        mark       = (x2_kv == x2_kv(index));
        [~,Vdiffh] = vext(Pts.y1_kv(mark),Pts.y2_kv(mark)); 
        
        
        if(order == 2)
            diffM      = Diff.DDy1; VdiffD     = Vdiffh.ddy1;
        elseif(order == 1)
            diffM      = Diff.Dy1; VdiffD     = Vdiffh.dy1;
        else
            diffM = eye(N1*N2); VdiffD = V(mark);
        end
        
        plot(Interp_r.pts1,Interp_r.InterPol*(diffM*V),'-'); hold on;
        plot(Pts.y1_kv(mark),diffM(mark,:)*V,'o'); hold on;
        plot(Pts.y1_kv(mark),VdiffD,'x');
        
        h = Pts.y2_kv(mark);
        title(['at theta = ',num2str(h(1))]);
    end

end