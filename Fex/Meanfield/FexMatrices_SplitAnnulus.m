function convStruct = FexMatrices_SplitAnnulus(optsPhys,IDC)
    
    global r_cutoff
    r_cutoff = optsPhys.V2.r_cutoff;

    params = optsPhys.V2;
    
    fstr      = (params.V2DV2);
    
    params     = rmfield(params,'V2DV2');
    paramNames = fieldnames(params);
    nParams    = length(paramNames);
    
    shapeAnn.RMin = params.LJsigma;
    shapeAnn.RMax = r_cutoff;
    if(isfield(optsPhys,'FexNum'))        
        shapeAnn.N    = optsPhys.FexNum.N;
    else            
        shapeAnn.N    = optsPhys.optsNum.PhysArea.Conv.N;
    end
    
    annulusArea       = Annulus(shapeAnn);

    shapeDisc.R = params.LJsigma;
    shapeDisc.N = shapeAnn.N;
    shapeDisc.volume = false;
    %diskArea    = Disc(shapeDisc);
    diskArea    = Sphere(shapeDisc);

    conv1 = IDC.ComputeConvolutionFiniteSupport(annulusArea,{fstr},IDC.Pts);
    conv2 = IDC.ComputeConvolutionFiniteSupport(diskArea,{fstr},IDC.Pts);

    convStruct.Conv =  conv1(:,:,2) + conv2(:,:,2);
    
    if(isfield(params,'epsilon'))
        convStruct.Conv = convStruct.Conv*params.epsilon;
    end

end