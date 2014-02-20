function SimulationHalfInfiniteCapillary

    disp('** Simulation Half Infinite Capillary **');

    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    N1 = 40; L1    = 1;
    N2 = 20; L2    = 1;
    vext  = @Vext7;
    
    x1PlotMax = Phys_to_Comp_1(3);
    y1Plot = Comp_to_Phys_1((-1:0.03:x1PlotMax)');
    y2Plot = Comp_to_Phys_2((-1:0.02:1)');    

    Maps = struct('PhysSpace1',@Comp_to_Phys_1,...
                  'PhysSpace2',@Comp_to_Phys_2,...
                  'CompSpace1',@Phys_to_Comp_1,...
                  'CompSpace2',@Phys_to_Comp_2);
    
     %Check Spectral/Spectral map in 2D:    
    [Pts,Diff,Int,Ind,Interp,Conv] = SpectralSpectral(Maps,N1,N2,y1Plot,y2Plot,@f2);                
    intBound = struct('y1_l',Comp_to_Phys_1(-1),...
                      'y1_u',Comp_to_Phys_1(1),...
                      'y2_l',Comp_to_Phys_2(-1),...
                      'y2_u',Comp_to_Phys_2(1));        
                       

    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound,'cart');    
    [VP]             = vext(Interp.pts1,Interp.pts2);           
                   
    %Check Differentiation
    vplot     = Interp.InterPol*V;        
    displayErrors(vplot,VP,V,Vdiff,Diff,'cart');
    
    %Check Interpolation
    doPlots_IP(Interp,V,VP);                    
    
    %Check Integration    
    display([' Error in Integration: ', num2str(Int(Int < inf)*V(Int < inf)-VInt)]);                
        
    %Show Convolution (no check)
    figure    
    doPlots_IP(Interp,Conv(:,Pts.y1_kv < inf)*V(Pts.y1_kv < inf));        
                                                                     
    %***************************************************************
    %   Mapping functions:
    %***************************************************************         
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys_1(x)
        [z,dz,dx,ddx,dddx,ddddx] = QuotientMap(x,L1,0);
    end
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys_2(xT)    
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xT,-L2,L2);
    end
    function xf = Phys_to_Comp_1(z)
         xf  = InvQuotientMap(z,L1,0);%  PosRay(z,L1);
    end
    function xf = Phys_to_Comp_2(z)           
        xf = z/L2;
    end
    %***************************************************************
    %   Auxiliary functions:
    %***************************************************************         
    function z = f2(x,y)
        z = exp(-(x.^2+y.^2));
    end
end