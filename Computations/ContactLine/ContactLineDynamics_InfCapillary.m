function [EX,res] = ContactLineDynamics_InfCapillary(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','InfCapillary_FMT','N',[40;40],... %
                       'L1',4,'y2Min',0,'y2Max',10,'N2bound',10);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',80,...
                       'y2Min',Phys_Area.y2Min,...
                       'y2Max',Phys_Area.y2Max,'N2',80);

    Sub_Area = struct('shape','Box','y1Min',-1,'y1Max',1,'N',[20,20],...
                      'y2Min',0,'y2Max',4);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);
	
    %V2Num    = struct('Fex','SplitAnnulus','N',[80,80]);
    %V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     
    
    plotTimes = struct('t_int',[0,5],'t_n',50);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'SubArea',Sub_Area,...
                     'FexNum',Fex_Num,...         
                     'plotTimes',plotTimes);%,...
                     %'V2Num',V2Num);
        
    V1 = struct('V1DV1','Vext_BarkerHenderson_InfCapillary',...
                'y2Min',Plot_Area.y2Min,'y2Max',Plot_Area.y2Max,...
                'epsilon_w',0.865*[1,1]);             
    
    optsPhys = struct('V1',V1,...%'V2',V2,...    
                      'kBT',0.75,...
                      'Dmu',0.0,...%'nParticlesS',5,...%
                      'nSpecies',1,'sigmaS',1); 

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    EX   = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    CL.ComputeEquilibrium(struct('solver','Newton'));    
%    EX.ComputeEquilibrium();%EX.optsPhys.rhoGas_sat);
    
    EX.IDC.plot(EX.GetRhoEq,'SC');
    EX.ComputeDynamics();
    res.fig_handles = EX.PlotDynamics();
end                 
