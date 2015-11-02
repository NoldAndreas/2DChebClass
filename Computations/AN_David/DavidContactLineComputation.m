function DavidContactLineComputation()

    %CL_LDA();
    CL_FMT();
    
    function CL_LDA()
        AddPaths();
        Phys_Area = struct('shape','HalfSpace_FMT','N',[30,30],...
                           'L1',4,'L2',5,...   
                           'y2wall',0.,...
                           'alpha_deg',45,...
                           'h',1);
                           %'N2bound',24,
                           %%'L2_AD',5.,...

        Plot_Area = struct('y1Min',-4,'y1Max',15,'N1',100,'N2',100,...
                           'y2Min',0.5,'y2Max',12.5);

        Sub_Area = struct('shape','Box','y1Min',-1,'y1Max',1,'N',[20,20],...
                          'y2Min',0,'y2Max',2);

        V2Num   = struct('Fex','SplitDisk','L',2,'L2',1.,'N',[20,20]);

        optsNum = struct('PhysArea',Phys_Area,...
                         'PlotAreaCart',Plot_Area,'SubArea',Sub_Area,...
                         'plotTimes',0:0.1:6,...
                         'V2Num',V2Num,...
                         'maxComp_y2',20);                                           

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',2.4);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

        optsPhys = struct('V1',V1,'V2',V2,...
                          'HSBulk','CarnahanStarling',...
                          'kBT',0.7,...                      
                          'Dmu',0.0,...
                          'sigmaS',1);

        CL = ContactLineHS(v2struct(optsPhys,optsNum));    
        CL.Preprocess();
        CL.ComputeEquilibrium();
        figure('Position',[0 0 600 600],'color','white');  CL.PlotContourResults({'plain'}); SaveFigure('LDA45degrees_contour'); %CL.IDC.plot(CL.GetRhoEq,'contour');
        figure('Position',[0 0 600 600],'color','white');  CL.IDC.plot(CL.GetRhoEq); SaveFigure('LDA45degrees');
        %figure;  CL.PlotContourResults();
    end

    function CL_FMT()
        
        AddPaths('ThesisNanoscale');   
               
        alpha_deg = 45;
        epw       = 1.155;        
        bounds1   = [-5 15];
        
        opts = v2struct(alpha_deg,epw,bounds1);            
        Job_ComputeContactAngle(opts);        
    end
    
    function config = GetStandardConfig(opts)                             
        
        config = ThesisNanoscale_GetStandardConfig(opts.alpha_deg,opts.epw);
                
        config.optsNum.PlotAreaCart       = struct('y1Min',opts.bounds1(1),'y1Max',opts.bounds1(2),...
                                                   'y2Min',0.5,'y2Max',12.5,...
                                                   'zMax',4,...
                                                    'N1',100,'N2',100);
    end

    function [nameEq,res] = Job_ComputeContactAngle(opts)
        
        %output:
        % - thetaY                  = Young contact angle
        % - int_DPI                 = -int_h0^\infty DP_I(h) dh
        % - error_int_DPI           = error of int_DPI
        % - thetaY_I                = contact angle obtained from integration of adsorption isotherm        
        % - error\Delta(\thetaY_I)  = error of thetaY_I
        % - int_DPII                = -int_{-\infty}^\infty DP_I(x) dx
        % - error_int_DPII          = error of int_DPII
        % - thetaY_II               = contact angle obtained from integration of DPII 
        % - error\Delta(\thetaY_II) = error of thetaY_II
        

        config = GetStandardConfig(opts);
        close all;

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium(struct('solver','Picard'));
               
        figure('Position',[0 0 600 600],'color','white');  CLT.PlotContourResults({'plain'}); SaveFigure('FMT45degrees_contour'); %CL.IDC.plot(CL.GetRhoEq,'contour');
        figure('Position',[0 0 600 600],'color','white');  CLT.IDC.plot(CLT.GetRhoEq); SaveFigure('FMT45degrees');
    end
       
end