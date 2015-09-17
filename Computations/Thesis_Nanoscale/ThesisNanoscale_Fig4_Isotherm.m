function ThesisNanoscale_Fig4_Isotherm()

    AddPaths('ThesisNanoscale');

	config = ThesisNanoscale_GetStandardConfig();
    config.optsNum.PhysArea.alpha_deg = 90;
    config.optsNum.PhysArea.N         = [1 250];
    config.optsNum.maxComp_y2         = -1;
    %config.optsNum.y1Shift           = 0;
    config.optsPhys.Dmu               = 0.03;
	    
    ComputeIsotherm(config);
    
    
    %For 45 degrees CA
    config.optsPhys.V1.epsilon_w = 1.155;
    ComputeIsotherm(config);
    
    %For 135 degrees CA
    config.optsPhys.V1.epsilon_w = 0.453;
    ComputeIsotherm(config);
    
    
    function ComputeIsotherm(config)
        CL = ContactLineHS(config);
        CL.Preprocess();    
        CL.ComputeAdsorptionIsotherm(600,'drying');    %wetting    
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(config.optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.2853);
        end
    end
 
end