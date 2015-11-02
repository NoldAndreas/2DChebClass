function data = SimulationInfCapillary_M1(N1,N2,L1,L2,vext)

    disp('** Simulation Inf Capillary M1 **');
    AddPaths();
    
    close all;
    
    %Initialization
    if(nargin == 0)
        N1 =  50;   N2 = 50; N = [N1;N2];
        L1 =  1;  y2Min = 0; y2Max = 2;
        vext  = @Vext_IN;
    end            
    
    IC     = InfCapillary_Interface(v2struct(L1,y2Min,y2Max,N));
    [Pts,Diff,Int,Ind] = IC.ComputeAll();
    Interp             = IC.ComputeInterpolationMatrix((-0.9:0.03:0.9)',...
                                (-1:0.02:1)',true,true);    
    
     %Check Spectral/Spectral map in 2D:         
    [V,Vdiff]   = vext(Pts.y1_kv,Pts.y2_kv);    
    [VP]        = vext(Interp.pts1,Interp.pts2);                          
    
    IC.plot(V);    
    
    %Check Interpolation        
    data.InterPol = max(abs(Interp.InterPol*V - VP));
    display([' Error in Interpolation: ', num2str(data.InterPol)]);
    
    %Check Differentiation
    vplot   = Interp.InterPol*V;        
    data    = displayErrorsPos(Pts,vplot,VP,V,Vdiff,Diff,'cart');
    
    subplot(2,2,1); IC.plot(Vdiff.dy2,'SC');
    subplot(2,2,2); IC.plot(Diff.Dy2*V,'SC');
    subplot(2,2,3); IC.plot(Diff.Dy2*V-Vdiff.dy2,'SC');
    
    figure;
    
    subplot(2,2,1); IC.plot(Vdiff.ddy2);
    subplot(2,2,2); IC.plot(Diff.DDy2*V,'SC');
    subplot(2,2,3); IC.plot(Diff.DDy2*V-Vdiff.ddy2,'SC');
       
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

    IC.plot(V,'SC');     title('Interpolation');        
    %***************************************************************
    %   Mapping functions:
    %***************************************************************         
    

    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function [h2,dh2,ddh2] = hP(y)
        y     = y - (y2Min+y2Max)/2;
        c     = 2.;
        h2    = c*y.^2;
        dh2   = 2*c*y;
        ddh2  = 2*c*ones(size(y));
        %h2   = c*y;
        %dh2  = c*ones(size(y));
        %ddh2 = zeros(size(y));
    end    
    function [V,VDiff] = Vext_IN(y1,y2)
        c2 = 0.1;
        
        [h2,dh2,ddh2] = hP(y2);

        z           = (y1-h2)/c2;
        V           = tanh(z);
        dVdy1       = (1-(tanh(z)).^2)/c2;
        dVdy2       = -(1-(tanh(z)).^2).*dh2/c2;
        dVddy1      = -2*tanh(z).*(1-(tanh(z).^2))/c2^2;
        
        %dVddy2      = -2*V.*(1-V.^2).*(dh2).^2/c2^2 - (1-V.^2).*ddh2/c2;
        
        dVddy2      = -(1-(tanh(z)).^2).*(...
                               2*tanh(z).*(dh2.^2)/c2^2+(ddh2)/c2);
        dVdy1dy2    = 2*tanh(z).*(1-(tanh(z)).^2).*dh2/c2^2;

         VDiff      = struct('dy1',dVdy1,'dy2',dVdy2,'ddy1',dVddy1,'ddy2',dVddy2,'dy1dy2',dVdy1dy2);

     end
   
end