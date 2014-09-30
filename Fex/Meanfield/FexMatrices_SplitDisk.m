function convStruct = FexMatrices_SplitDisk(optsPhys,IDC)
  
    params = optsPhys.V2;
    
    fstr      = (params.V2DV2);
    
    params     = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams    = length(paramNames);
    
    shapeAnn.RMin = params.LJsigma;
    if(isfield(optsPhys,'FexNum'))
        shapeAnn.L    = optsPhys.FexNum.L;
        shapeAnn.N    = optsPhys.FexNum.N;
    else    
        shapeAnn.L    = optsPhys.optsNum.PhysArea.Conv.L;
        shapeAnn.N    = optsPhys.optsNum.PhysArea.Conv.N;
    end
    annulusArea       = InfAnnulus(shapeAnn);

    shapeDisc.R = params.LJsigma;
    shapeDisc.N = shapeAnn.N;
    diskArea    = Disc(shapeDisc);

    conv1 = IDC.ComputeConvolutionFiniteSupport(annulusArea,{fstr},IDC.Pts);
    conv2 = IDC.ComputeConvolutionFiniteSupport(diskArea,{fstr},IDC.Pts);

    convStruct.Conv = conv1(:,:,2) + conv2(:,:,2);
    
    if(isfield(params,'epsilon'))
        convStruct.Conv = convStruct.Conv*params.epsilon;
    end

end