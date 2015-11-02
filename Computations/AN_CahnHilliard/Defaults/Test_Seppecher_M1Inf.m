function [optsNum,optsPhys] = Test_Seppecher_M1Inf()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'N1',50,...
                       'y2Min',0,'y2Max',20,'N2',20,'L1',12); %y1Max = 20
    
    Plot_Area = struct('y1Min',-20,'y1Max',20,'N1',120,...
                       'y2Min',0,'y2Max',Phys_Area.y2Max,'N2',40,...
                       'N1Vecs',40,'N2Vecs',6,'Hy2',3);
                   
    Sub_Area  = struct('y1Min',-5,'y1Max',5,'N1',50,...
                       'y2Min',0,'y2Max',Phys_Area.y2Max,'N2',30);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'DDFTCode','Seppecher_M1Inf',...
                     'plotTimes',0:0.1:5);                            
                      
    optsPhys = struct('theta',pi/2,'g',0,...
                      'D_A',0.,'rho_m',2,'nu',10,'Ca',0.02,... 
                     'nParticles',0,'UWall',-1);
                 
    set(0,'defaultaxesfontsize',20);
    set(0,'defaultlinelinewidth',1.);
    set(0,'defaulttextfontsize',15.);
    %set(0,'defaultaxeswidth',1);
                 
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);
    
    close all;
end                 

