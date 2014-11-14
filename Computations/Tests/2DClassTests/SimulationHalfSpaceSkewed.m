function data = SimulationHalfSpaceSkewed(N1,N2,L1,L2)

    disp('** Simulation Half Space **');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    if(nargin == 0)        
        VextC  = @Vext15;        
        Phys_Area = struct('y1Min',-inf,'y1Max',inf,...
                           'L1',3,'L2',2,'N',[20;20],...
                           'y2Min',0,'alpha',pi/4);

        PlotArea       = struct('y1Min',-5,'y1Max',5,'L1',3,'L2',2,...
                           'N2',100,'N1',100,'y2Min',0,'y2Max',4);    

        Phys_Area.Conv      = struct('L',3,'L2',1,'N',[20,20]);%'ep2Conv',0.1                       
        
    else
        N = [N1;N2];
    end        
    
    HS                        = HalfSpaceSkewed(Phys_Area);%v2struct(L1,L2,N));
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(PlotArea);        
    
    %[Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(PlotArea);        
    
    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv);
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
    Conv = HS.ComputeConvolutionMatrix(@f2,true);
    yPtsCheck          = [ 20 10; 0 2 ; 2 0 ; -10 0];
    HS.TestConvolutionMatrix(yPtsCheck,@f2);
        
    %Check Convolution at infinity
    [~,intCheck] = f2(0);
    int          = Conv*ones(size(V));
    disp(['Error of convolution at infinity: ',num2str(max(abs(int(Pts.y2_kv==inf)-intCheck)))]);
    
    fP = V;%f1(Pts.y1_kv,Pts.y2_kv);
    %data.Conv = max(abs(Interp.InterPol*(Conv*fP) - fConv(Interp.pts1,Interp.pts2)));
    %display([' Error in Convolution: ', num2str(data.Conv)]);
    
%    data.N1 = N(1); data.N2 = N(2);
%    data.Interp = Interp; data.f = V;    
    
    %******** Plotting **********
    figure
    set(gcf,'Color','white'); %Set background color    
    
    subplot(1,2,1);
    HS.plot(V,'SC');    
    title('Interpolation');    
    pbaspect([1 1 1]);
    
    subplot(1,2,2);
    HS.plot(Conv*fP,'SC');    
    title('Convolution');
    pbaspect([1 1 1]);

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [z,int] = f2(r)        
        z       = exp(-r.^2);
        int     = pi;
    end

    function [V,Vdiff,VInt] = vext(y1t,y2t)     
        ptsCart = HS.GetCartPts(y1t,y2t);
        %[y1,y2]        = SkewedGrid(y1t,y2t,Phys_Area.alpha);        
        [V,Vdiff,VInt] = VextC(ptsCart.y1_kv,ptsCart.y2_kv);
    end
end