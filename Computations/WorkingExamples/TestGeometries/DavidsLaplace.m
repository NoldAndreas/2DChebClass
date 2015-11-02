function data = DavidsLaplace()
%************************************************
%  define mu_s = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') + V_ext - mu
% Equilibrium:
% (EQ 1)  mu_s      = 0
% (EQ 2)  int(rho)  = noParticles
% Dynamics:
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC 1)  n*grad(rho) = 0
%************************************************


    Phys_Area = struct('y1Min',0,'y1Max',10,'N1',25,...
                       'y2Min',0,'y2Max',10,'N2',25);
	Plot_Area = struct('y1Min',1,'y1Max',10,'N1',100,...
                       'y2Min',0.1,'y2Max',9.9,'N2',100);
    
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',1,'y1Max',3,'N1',20,...
                      'y2Min',1,'y2Max',2,'N2',20);
    
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:3,...
                     'DDFTCode','DDFT_DiffusionBox');    
                      
    disp(optsNum.DDFTCode);
 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;    
    [N1,N2,PhysArea,~,y1Plot,y2Plot,~] = LoadNumData(optsNum);    

    Maps = struct('PhysSpace1',@Comp_to_Phys1,...
                  'PhysSpace2',@Comp_to_Phys2,...
                  'CompSpace1',@Phys_to_Comp1,...
                  'CompSpace2',@Phys_to_Comp2);                
              
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    
    [~,Diff,~,Ind,Interp] = SpectralSpectral(Maps,N1,N2,y1Plot,y2Plot);              

    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************        
    x_ic    = fsolve(@f,zeros(N1*N2,1));
    doPlots_IP(Interp,x_ic);
    
    SaveToFile('DavidsPicData',v2struct(x_ic,Interp),getResultsPath());

    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function y = f(x)
        %solves for T*log*rho + Vext                
        y            = Diff.Lap*x;
        y(Ind.bound) = x(Ind.bound);
        %y(Ind.bound) = Ind.normal*(Diff.grad*x);
        y(Ind.right) = x(Ind.right) -1;
        y(Ind.left)  = x(Ind.left);
    end
        
    %***************************************************************
    %Mapping functions:
    %***************************************************************
    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys1(xf)
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xf,PhysArea.y1Min,PhysArea.y1Max);
    end
    function xf = Phys_to_Comp1(z)
        xf = InvLinearMap(z,PhysArea.y1Min,PhysArea.y1Max);        
    end

    function [z,dz,dx,ddx,dddx,ddddx] = Comp_to_Phys2(xT)    
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xT,PhysArea.y2Min,PhysArea.y2Max);
    end

    function xf = Phys_to_Comp2(z)           
        xf = InvLinearMap(z,PhysArea.y2Min,PhysArea.y2Max);        
    end    

end