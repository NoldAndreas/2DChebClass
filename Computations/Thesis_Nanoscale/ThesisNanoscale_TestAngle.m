function ThesisNanoscale_TestAngle()

    config = ThesisNanoscale_GetStandardConfig(30,0.856);
    config.optsNum.PhysArea.N         = [1 20];
    config.optsNum.maxComp_y2         = -1;
    %config.optsNum.y1Shift           = 0;
    config.optsPhys.Dmu               = 0.0;
    
    config.optsNum.V2Num.N = [20,20];
    config.optsNum.V2Num.N1disc = 20;
    config.optsNum.V2Num.N2disc = 20;
    
    CL = ContactLineHS(config);
    CL.Preprocess();    
	CL.ComputeAdsorptionIsotherm(600,'drying');    %wetting    
    CL.FittingAdsorptionIsotherm([10 14],1)
    if(config.optsPhys.kBT == 0.75)
        CL.SumRule_AdsorptionIsotherm(0.2853);
    end

end