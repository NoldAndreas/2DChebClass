function SimulationBoxFD

    disp('** SimulationBoxFD **');
    if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %Initialization
    L1    = 1;
    L2    = 1;
    vext  = @Vext5;
    
    rPlot = Comp_to_Phys_R((-1:0.03:1)');
    tPlot = Comp_to_Phys_T((-1:0.02:1)');    

    Maps = struct('PhysSpace1',@Comp_to_Phys_R,...
                  'PhysSpace2',@Comp_to_Phys_T,...
                  'CompSpace1',@Phys_to_Comp_R,...
                  'CompSpace2',@Phys_to_Comp_T);
         
    [Pts,Diff,Int,Ind,Interp] = FDFD(Maps,50,50,rPlot,tPlot);                
    intBound = struct('y1_l',Comp_to_Phys_R(-1),...
                      'y1_u',Comp_to_Phys_R(1),...
                      'y2_l',Comp_to_Phys_T(-1),...
                      'y2_u',Comp_to_Phys_T(1));        
                       

    [V,Vdiff,VInt]   = vext(Pts.y1_kv,Pts.y2_kv,intBound,'cart');    
    [VP]             = vext(Interp.pts1,Interp.pts2);           
                   
    %Check Interpolation    
    doPlots_IP(Interp,V,VP);                    
    
    %Check Differentiation and Interpolation    
    vplot     = Interp.InterPol*V;        
    displayErrors(vplot,VP,V,Vdiff,Diff,'cart');    
    
    %Check Integration
    display([' Error in Integration: ', num2str(Int*V-VInt)]);                
                                            
    
    %***************************************************************
    %   Mapping functions:
    %***************************************************************         
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys_R(xR)
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xR,-L1,L1);
    end
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys_T(xT)    
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xT,-L2,L2);
    end
    function xf = Phys_to_Comp_R(z)
         xf  = z/L1;
    end
    function xf = Phys_to_Comp_T(z)           
        xf = z/L2;
    end
function doPlots_IP(Interp,V,VP)        
    
    nSpecies=size(V,2);
    if(nSpecies == 1)
        nCol = 1;
    else
        nCol = 2;
    end    
    nRows=ceil(nSpecies/nCol);
    
    y1 = Interp.pts1;
    y2 = Interp.pts2;   
    
    if(length(V) == length(Interp.pts1))
        z = V;
    else
        z = real(Interp.InterPol*V);        
    end

    y1M     = reshape(y1,Interp.Nplot2,Interp.Nplot1);
    y2M     = reshape(y2,Interp.Nplot2,Interp.Nplot1);
    
    xl = [(min(y1)-0.5) (max(y1)+0.5)];
    yl = [(min(y2)-0.5) (max(y2)+0.5)];
        
    if(nargin == 2)
        
        if(nSpecies >= 2)
        
            for iSpecies=1:nSpecies
                subplot(nRows,nCol,iSpecies);
                mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));         
                
                SetXYLabel();                                
                title(['Species ' num2str(iSpecies)]);
            end
        else
            mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));         
            SetXYLabel();
        end

    elseif(nargin == 3)
        
        if(length(VP) == length(Interp.pts1))
            zP = VP;
        else
            zP = real(Interp.InterPol*VP);        
        end    
        
        subplot(2,1,1)        
        mesh(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1));        
        SetXYLabel();
            
        subplot(2,1,2)        
        mesh(y1M,y2M,reshape(z-zP,Interp.Nplot2,Interp.Nplot1));                           
        SetXYLabel();
    end
    pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1/2*min((xl(2)-xl(1)),(yl(2)-yl(1)))]);    
	set(gca,'fontsize',15);        
    
    function SetXYLabel()
        h = xlabel('$y_1$');
        set(h,'Interpreter','Latex'); 
        set(h,'fontsize',20);
        
        h = ylabel('$y_2$');
        set(h,'Interpreter','Latex'); 
        set(h,'fontsize',20);
    end
        
    
end


end