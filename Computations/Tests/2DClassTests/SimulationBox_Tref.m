function data = SimulationBox_Tref(N1,N2,L1,L2,vext)

    disp('** Simulation Box Tref **');
    AddPaths();    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  30;   N2 = 30;
        PlotArea = struct('y1Min',-1,'y1Max',1,'N1',70,...
                          'y2Min',-1,'y2Max',1,'N2',70);
                
        vext  = @Vext16;
    end
    
    TS                        = BoxTrefSpectralSpectral(N1,N2);
    [Pts,Diff,Int,Ind,Interp] = TS.ComputeAll(PlotArea);
            
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP,VPDiff]        = vext(Interp.pts1,Interp.pts2);  
    
    subplot(1,2,1); TS.plot(V,'SC');    
    
    [V,Pts,Diff,Int,Ind,Interp] = TS.UpdatePadeValues(V,PlotArea);
                                            
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);                          
        
    subplot(1,2,1); TS.plot(Vdiff.ddy2,'SC');
    subplot(1,2,2); TS.plot(Vdiff.ddy1,'SC');
    
    subplot(1,2,2); 
    TS.plot(V,'SC'); 
    TS.PlotLineOfPoles(V);
     
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    figure;
    subplot(2,2,1); TS.plot(Vdiff.ddy2,'SC');
    subplot(2,2,2); TS.plot(Diff.DDy2*V,'SC');
    subplot(2,2,3); TS.plot(Diff.DDy2*V-Vdiff.ddy2,'SC');
    
    figure;
    subplot(2,2,1); TS.plot(Vdiff.dy1dy2,'SC');
    subplot(2,2,2); TS.plot(Diff.Dy1Dy2*V,'SC');
    subplot(2,2,3); TS.plot(Diff.Dy1Dy2*V-Vdiff.dy1dy2,'SC');
        
    %Check Integration
    %data.Int = abs(Int*V-VInt);
    %display([' Error in Integration: ', num2str(data.Int)]);                
        
    %Check Convolution
    %fP = f1(Pts.y1_kv,Pts.y2_kv);
    %data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    %display([' Error in Convolution: ', num2str(data.Conv)]);
    
    data.N1 = N1; data.N2 = N2;
%    data.Interp = Interp; data.f = V;    
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    TS.plot(V);    
    title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);
    TS.plot(fConv(Interp.pts1,Interp.pts2));
    title('Convolution');
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
end